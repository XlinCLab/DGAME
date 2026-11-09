"""PyQt6 window for setting up a DGAME experiment config."""

import codecs
import os
import re
import signal
import sys
import tempfile
from pathlib import Path

import yaml
from PyQt6.QtCore import QProcess, QProcessEnvironment, QSize, Qt, QTimer
from PyQt6.QtGui import QFontDatabase, QKeySequence, QShortcut, QTextCursor
from PyQt6.QtWidgets import (QApplication, QCheckBox, QComboBox, QFileDialog,
                             QFrame, QHBoxLayout, QLabel, QLineEdit,
                             QMainWindow, QMessageBox, QPlainTextEdit,
                             QPushButton, QScrollArea, QTabWidget, QTreeWidget,
                             QTreeWidgetItem, QVBoxLayout, QWidget)

from dgame.pipeline import JULIA_STEPS, R_STEPS
from experiment.gui.config_model import (BOOL, STR_LIST, VERSION_KEY, YAML,
                                         ConfigModel, Field,
                                         available_versions, dependency_status,
                                         discover_subjects, load_defaults,
                                         value_to_text)

# Fields shown on the "Setup" tab (the rest are on the "Advanced" tab)
SETUP_KEYS = [
    "experiment.id",
    "experiment.input_data.root",
    "experiment.outdir",
    "experiment.subjects",
    "experiment.objects",
    "experiment.fillers",
]
SETUP_LABELS = {
    "experiment.id": "Experiment ID",
    "experiment.input_data.root": "Input data folder",
    "experiment.outdir": "Output folder",
    "experiment.subjects": "Subjects",
    "experiment.objects": "Target words",
    "experiment.fillers": "Filler words",
}
SETUP_HINTS = {
    "experiment.id": "Name of the output subfolder; a timestamp is used if empty.",
    "experiment.input_data.root": "Folder containing the preproc/ and recordings/ input data.",
    "experiment.outdir": "Results are written to a subfolder named by the experiment ID. Defaults to ./out.",
    "experiment.subjects": "One subject ID per line. Leave empty to process all subjects found in the input data.",
    "experiment.objects": "One word per line, written as the nouns' lemmas (e.g. Kerze).",
    "experiment.fillers": "One word per line, written as the nouns' lemmas.",
}
ERROR_STYLE = "border: 1px solid #d33;"
RUN_SCRIPT = Path(__file__).resolve().parents[2] / "run_dgame_analysis.py"
ANSI_ESCAPE_REGEX = re.compile(r"\x1b\[[0-9;?]*[ -/]*[@-~]")
STOP_TIMEOUT_MS = 5000
LOG_FILE_MESSAGE_REGEX = re.compile(r"Experiment run log file: (.+)")
CRASH_REPORT_MARKER = "Fatal Python error"
MAX_CRASH_REPORT_CHARS = 10000


def describe_signal(number: int) -> str:
    """Description of the signal which killed a process, e.g. "SIGSEGV: Segmentation fault: 11"."""
    try:
        return f"{signal.Signals(number).name}: {signal.strsignal(number)}"
    except ValueError:
        return f"signal {number}"


class ConfigWindow(QMainWindow):
    def __init__(self, model: ConfigModel, required_fields: list[str] | None = None):
        super().__init__()
        self.setWindowTitle("DGAME experiment setup")
        self.resize(1100, 800)
        self.model = model
        self.required_fields = required_fields or []
        self.process: QProcess | None = None  # the running analysis, if any
        self.config_file: str | None = None  # temporary config file of the running analysis
        self.run_message = ""  # outcome of the last analysis run
        self.stopping = False

        self.tabs = QTabWidget()
        self.status = QLabel()
        self.status.setWordWrap(True)
        self.log_view = QPlainTextEdit()
        self.log_view.setReadOnly(True)
        self.log_view.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.log_view.setMaximumBlockCount(50000)
        self.log_view.setFont(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont))
        self.log_view.setPlaceholderText("The log of the analysis will appear here once it is running.")
        self.run_button = QPushButton("Run analysis")
        self.run_button.clicked.connect(self.toggle_analysis)
        self.quit_button = QPushButton("Quit")
        self.quit_button.clicked.connect(self.close)
        QShortcut(QKeySequence.StandardKey.Quit, self, self.close)
        self.load_button, self.save_button = QPushButton("Load config…"), QPushButton("Save config…")
        self.load_button.clicked.connect(self.load_from_file)
        self.save_button.clicked.connect(self.save_to_file)

        buttons = QHBoxLayout()
        buttons.addWidget(self.load_button)
        buttons.addWidget(self.save_button)
        buttons.addWidget(self.status, stretch=1)
        buttons.addWidget(self.run_button)
        buttons.addWidget(self.quit_button)
        central = QWidget()
        layout = QVBoxLayout(central)
        layout.addWidget(self.tabs)
        layout.addLayout(buttons)
        self.setCentralWidget(central)
        self.rebuild()

    # --- Building the tabs ---

    def rebuild(self) -> None:
        """(Re)create all tabs from the current model."""
        current = self.tabs.currentIndex()
        self.tabs.clear()
        self.tabs.addTab(self.build_setup_tab(), "Setup")
        self.tabs.addTab(self.build_pipeline_tab(), "Pipeline")
        self.tabs.addTab(self.build_advanced_tab(), "Advanced")
        self.tabs.addTab(self.log_view, "Log")
        self.tabs.setCurrentIndex(max(current, 0))
        self.update_status()

    def label_for(self, key: str) -> str:
        return SETUP_LABELS.get(key, key) + (" *" if self.is_required(key) else "")

    def is_required(self, key: str) -> bool:
        return any(re.search(p, key) for p in self.required_fields)

    def build_setup_tab(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setSpacing(4)
        versions = available_versions()
        version_box = QComboBox()
        version_box.addItems([f"DGAME {v}" for v in versions])
        current = str(self.model.values[VERSION_KEY])
        if current in versions:
            version_box.setCurrentIndex(versions.index(current))
        version_box.currentIndexChanged.connect(lambda i: self.change_version(versions[i]))
        version_row = QHBoxLayout()
        version_row.addWidget(QLabel("Experiment version"))
        version_row.addWidget(version_box)
        version_row.addStretch()
        layout.addLayout(version_row)
        layout.addSpacing(12)

        # NB: each field is a label above its editor and hint (rather than a form layout row), as wrapped
        # hint labels are only given enough height when they are laid out directly in a vertical layout
        for key in SETUP_KEYS:
            if key not in self.model.fields:
                continue
            label = QLabel(self.label_for(key))
            font = label.font()
            font.setBold(True)
            label.setFont(font)
            layout.addWidget(label)
            editor = self.make_editor(self.model.fields[key], browse=self.model.fields[key].is_dir)
            layout.addWidget(editor)
            if key == "experiment.subjects":
                find = QPushButton("Find subjects in input folder")
                find.clicked.connect(lambda _=False, e=editor: self.find_subjects(e))
                layout.addWidget(find, alignment=Qt.AlignmentFlag.AlignLeft)
            hint = QLabel(SETUP_HINTS.get(key, ""))
            hint.setStyleSheet("color: gray;")
            hint.setWordWrap(True)
            layout.addWidget(hint)
            layout.addSpacing(12)
        required_note = QLabel("* required")
        required_note.setStyleSheet("color: gray;")
        layout.addWidget(required_note)
        layout.addStretch()
        # Scrollable, so that rows keep the height their font needs instead of being squeezed in a small window
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setWidget(page)
        return scroll

    def build_pipeline_tab(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        tree = QTreeWidget()
        tree.setColumnCount(3)
        tree.setHeaderLabels(["Analysis step", "Requires", "Description"])
        tree.setRootIsDecorated(False)
        available = dependency_status()
        for f in self.model.fields.values():
            if not f.is_step_switch:
                continue
            step_id = f.path[2]
            needs = "Julia" if step_id in JULIA_STEPS else "R" if step_id in R_STEPS else ""
            item = QTreeWidgetItem(["", needs + ("" if not needs or available[needs] else " (not found)"), self.model.comments.get(".".join(f.path[:3]), "")])
            if needs and not available[needs]:
                item.setForeground(1, Qt.GlobalColor.red)
            tree.addTopLevelItem(item)
            # A checkbox widget per row (rather than item check states, which the native macOS style only draws for the first row)
            box = QCheckBox(step_id)
            box.setChecked(bool(self.model.values[f.key]))
            box.toggled.connect(lambda checked, key=f.key: self.edited(key, None, checked))
            tree.setItemWidget(item, 0, box)
        for column in (0, 1):
            tree.resizeColumnToContents(column)
        layout.addWidget(tree)
        row = QHBoxLayout()
        for text, state in (("Enable all", True), ("Disable all", False)):
            button = QPushButton(text)
            button.clicked.connect(lambda _=False, s=state: self.set_all_steps(tree, s))
            row.addWidget(button)
        row.addStretch()
        layout.addLayout(row)
        note = QLabel("Steps depend on the outputs of earlier steps; disabling a step only works if its outputs already exist.")
        note.setStyleSheet("color: gray;")
        layout.addWidget(note)
        return page

    def build_advanced_tab(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        search = QLineEdit()
        search.setPlaceholderText("Filter parameters…")
        search.setClearButtonEnabled(True)
        tree = QTreeWidget()
        tree.setColumnCount(3)
        tree.setHeaderLabels(["Parameter", "Value", "Description"])
        groups: dict[tuple, QTreeWidgetItem] = {}

        def group_item(path: tuple) -> QTreeWidgetItem:
            if path not in groups:
                item = QTreeWidgetItem([path[-1]])
                (group_item(path[:-1]) if len(path) > 1 else tree.invisibleRootItem()).addChild(item)
                font = item.font(0)
                font.setBold(True)
                item.setFont(0, font)
                groups[path] = item
            return groups[path]

        for key, f in self.model.fields.items():
            if key in SETUP_KEYS or key == VERSION_KEY or f.is_step_switch:
                continue
            item = QTreeWidgetItem([f.label, "", f.help])
            item.setToolTip(0, key)
            item.setToolTip(2, f.help)
            item.setData(0, Qt.ItemDataRole.UserRole, key)
            group_item(f.path[:-1]).addChild(item)
            editor = self.make_editor(f, browse=f.is_dir, text_lines=2)
            tree.setItemWidget(item, 1, editor)
            multiline = isinstance(editor, QPlainTextEdit)
            item.setSizeHint(1, QSize(300, editor.maximumHeight() + 6) if multiline else editor.sizeHint())
        tree.setColumnWidth(0, 300)
        tree.setColumnWidth(1, 300)
        search.textChanged.connect(lambda text: self.filter_tree(tree, text))
        layout.addWidget(search)
        layout.addWidget(tree)
        return page

    # --- Editors ---

    def make_editor(self, f: Field, browse: bool = False, text_lines: int = 3) -> QWidget:
        key = f.key
        value = self.model.values[key]
        if f.kind == BOOL:
            box = QCheckBox()
            box.setChecked(bool(value))
            box.toggled.connect(lambda checked: self.edited(key, None, checked))
            return box
        text = value_to_text(f.kind, value)
        if f.kind in (STR_LIST, YAML):
            edit = QPlainTextEdit(text)
            edit.setFixedHeight(self.text_height(edit, text_lines))
            edit.textChanged.connect(lambda: self.edited(key, edit, edit.toPlainText()))
            return edit
        edit = QLineEdit(text)
        edit.setMinimumHeight(edit.fontMetrics().height() + 12)
        if f.kind == "optional_str":
            edit.setPlaceholderText("(not set)" if key != "experiment.outdir" else "out")
        edit.textChanged.connect(lambda t: self.edited(key, edit, t))
        if not browse:
            return edit
        button = QPushButton("Browse…")
        button.clicked.connect(lambda: self.browse_dir(edit))
        row = QWidget()
        row_layout = QHBoxLayout(row)
        row_layout.setContentsMargins(0, 0, 0, 0)
        row_layout.addWidget(edit, stretch=1)
        row_layout.addWidget(button)
        return row

    @staticmethod
    def text_height(edit: QPlainTextEdit, lines: int) -> int:
        """Height of a multi-line text editor showing a number of lines in its font."""
        margins = 2 * edit.frameWidth() + 2 * int(edit.document().documentMargin())
        return lines * edit.fontMetrics().lineSpacing() + margins + 8

    def edited(self, key: str, widget: QWidget | None, value) -> None:
        if widget is None:
            self.model.set_value(key, value)
        else:
            error = self.model.set_text(key, value)
            widget.setStyleSheet(ERROR_STYLE if error else "")
            widget.setToolTip(error or "")
        self.update_status()

    def browse_dir(self, edit: QLineEdit) -> None:
        directory = QFileDialog.getExistingDirectory(self, "Choose folder", edit.text() or os.getcwd())
        if directory:
            edit.setText(directory)

    # --- Actions ---

    def find_subjects(self, editor: QWidget) -> None:
        found = discover_subjects(self.model)
        if not found:
            QMessageBox.information(self, "No subjects found", "No subject folders were found in <input folder>/recordings/xdf.")
            return
        editor.setPlainText("\n".join(found))

    @staticmethod
    def set_all_steps(tree: QTreeWidget, enabled: bool) -> None:
        for i in range(tree.topLevelItemCount()):
            tree.itemWidget(tree.topLevelItem(i), 0).setChecked(enabled)

    @staticmethod
    def filter_tree(tree: QTreeWidget, text: str) -> None:
        text = text.strip().lower()

        def visit(item: QTreeWidgetItem) -> bool:
            if item.childCount() == 0:
                haystack = " ".join([item.text(0), item.text(2), str(item.data(0, Qt.ItemDataRole.UserRole) or "")]).lower()
                visible = text in haystack
            else:
                visible = any([visit(item.child(i)) for i in range(item.childCount())])
            item.setHidden(not visible)
            return visible

        for i in range(tree.topLevelItemCount()):
            visit(tree.topLevelItem(i))
        # Show matches while filtering; collapse everything again when the filter is cleared
        tree.expandAll() if text else tree.collapseAll()

    def change_version(self, version: str) -> None:
        if str(self.model.values[VERSION_KEY]) == version:
            return
        defaults, comments = load_defaults(version)
        self.model.rebase(defaults, comments)
        self.rebuild()

    def load_from_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Load config", os.getcwd(), "YAML files (*.yml *.yaml)")
        if not path:
            return
        try:
            with open(path) as f:
                config = yaml.safe_load(f) or {}
            version = str(config.get("experiment", {}).get("dgame_version", self.model.values[VERSION_KEY]))
            defaults, comments = load_defaults(version if version in available_versions() else self.model.values[VERSION_KEY])
            model = ConfigModel(defaults, comments)
            model.load(config)
        except Exception as exc:
            QMessageBox.critical(self, "Could not load config", str(exc))
            return
        self.model = model
        self.rebuild()

    def save_to_file(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save config", os.getcwd(), "YAML files (*.yml *.yaml)")
        if path:
            with open(path, "w") as f:
                yaml.safe_dump(self.model.to_dict(only_changes=True), f, sort_keys=False, allow_unicode=True)

    def update_status(self) -> None:
        if self.is_running():
            self.run_button.setEnabled(True)
            self.run_button.setToolTip("")
            self.status.setText("Analysis running…")
            return
        problems = self.model.problems(self.required_fields)
        self.run_button.setEnabled(not problems)
        self.run_button.setToolTip("\n".join(problems))
        if problems:
            self.status.setText("<span style='color:#c33'>Cannot run yet:</span> " + "; ".join(problems[:3]) + (f" (+{len(problems) - 3} more)" if len(problems) > 3 else ""))
        else:
            self.status.setText(self.run_message or "Ready to run.")

    # --- Running the analysis ---
    # The analysis runs in a separate process (as when run from the command line), so that it cannot freeze or crash
    # this window, and so that it can be stopped. Its log output is shown in the "Log" tab.

    def is_running(self) -> bool:
        return self.process is not None

    def analysis_command(self, config_path: str) -> tuple[str, list[str]]:
        """Program and arguments to run the analysis from a config file."""
        return sys.executable, ["-u", str(RUN_SCRIPT), "--config", config_path]

    def toggle_analysis(self) -> None:
        self.stop_analysis() if self.is_running() else self.start_analysis()

    def confirm_overwrite(self, config: dict) -> str | None:
        """Ask whether to overwrite the output of a previous run with the same experiment ID.
        Returns the answer to give the analysis' prompt ("y" or "n"; "" if there is no prompt), or None to cancel."""
        experiment = config.get("experiment", {})
        experiment_id = (experiment.get("id") or "").strip()
        if not experiment_id:
            return ""  # the ID defaults to a timestamp, so the output folder is new
        output_dir = os.path.join(os.path.abspath(experiment.get("outdir") or "out"), experiment_id)
        if not os.path.exists(output_dir):
            return ""
        reply = QMessageBox.question(
            self, "Output folder exists",
            f"The output folder\n{output_dir}\nalready exists. Overwrite the previous results?\n\n"
            "Choose No to keep them and write to a new folder with a timestamp.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.No,
        )
        return {QMessageBox.StandardButton.Yes: "y", QMessageBox.StandardButton.No: "n"}.get(reply)

    def start_analysis(self) -> None:
        config = self.model.to_dict()
        answer = self.confirm_overwrite(config)
        if answer is None:
            return
        handle, self.config_file = tempfile.mkstemp(prefix="dgame_gui_config_", suffix=".yml")
        with os.fdopen(handle, "w") as f:
            yaml.safe_dump(config, f, sort_keys=False, allow_unicode=True)

        process = QProcess(self)
        process.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        environment = QProcessEnvironment.systemEnvironment()
        environment.insert("PYTHONUNBUFFERED", "1")
        process.setProcessEnvironment(environment)
        process.readyReadStandardOutput.connect(self.read_output)
        process.finished.connect(self.analysis_finished)
        process.errorOccurred.connect(self.process_error)
        process.started.connect(lambda: self.send_answer(process, answer))
        self.process, self.stopping, self.run_message = process, False, ""
        self.decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")
        self.log_view.clear()
        self.set_running(True)
        self.tabs.setCurrentWidget(self.log_view)
        program, arguments = self.analysis_command(self.config_file)
        process.start(program, arguments)

    @staticmethod
    def send_answer(process: QProcess, answer: str) -> None:
        """Answer the analysis' prompt about overwriting the output folder (if it asks) and close its input."""
        if answer:
            process.write((answer + "\n").encode())
        process.closeWriteChannel()

    def stop_analysis(self) -> None:
        process = self.process
        if process is None:
            return
        self.stopping = True
        self.append_log("\n[Stopping the analysis…]\n")
        process.terminate()
        QTimer.singleShot(STOP_TIMEOUT_MS, lambda: process.kill() if process.state() != QProcess.ProcessState.NotRunning else None)

    def set_running(self, running: bool) -> None:
        for index in range(self.tabs.count()):
            if self.tabs.widget(index) is not self.log_view:
                self.tabs.setTabEnabled(index, not running)
        self.load_button.setEnabled(not running)
        self.save_button.setEnabled(not running)
        self.run_button.setText("Stop" if running else "Run analysis")
        self.update_status()

    def read_output(self) -> None:
        if self.process is not None:
            self.append_log(self.decoder.decode(bytes(self.process.readAllStandardOutput())))

    def append_log(self, text: str) -> None:
        """Add text to the log, treating a carriage return as overwriting the current line (as progress bars do)."""
        text = ANSI_ESCAPE_REGEX.sub("", text)
        scrollbar = self.log_view.verticalScrollBar()
        follow = scrollbar.value() >= scrollbar.maximum() - 2
        cursor = QTextCursor(self.log_view.document())
        cursor.movePosition(QTextCursor.MoveOperation.End)
        for part in re.split(r"(\r\n|\n|\r)", text):
            if part in ("\n", "\r\n"):
                cursor.insertBlock()
            elif part == "\r":
                cursor.movePosition(QTextCursor.MoveOperation.StartOfBlock, QTextCursor.MoveMode.KeepAnchor)
                cursor.removeSelectedText()
            elif part:
                cursor.insertText(part)
        if follow:
            scrollbar.setValue(scrollbar.maximum())

    def process_error(self, error: QProcess.ProcessError) -> None:
        if error == QProcess.ProcessError.FailedToStart:
            self.append_log("Could not start the analysis.\n")
            self.analysis_finished(-1, QProcess.ExitStatus.CrashExit)

    def analysis_finished(self, exit_code: int, exit_status: QProcess.ExitStatus) -> None:
        if self.process is None:
            return
        self.read_output()
        self.append_log(self.decoder.decode(b"", final=True))
        if self.stopping:
            self.run_message = "Analysis stopped."
        elif exit_status == QProcess.ExitStatus.CrashExit:
            # On Unix, the "exit code" of a crashed process is the number of the signal which killed it
            self.run_message = f"<span style='color:#c33'>The analysis crashed</span> ({describe_signal(exit_code)}; see the log)."
            self.show_crash_report()
        elif exit_code != 0:
            self.run_message = f"<span style='color:#c33'>The analysis failed</span> (exit code {exit_code}; see the log)."
        else:
            self.run_message = "Analysis finished successfully."
        self.append_log(f"\n[{re.sub('<[^>]+>', '', self.run_message)}]\n")
        self.process.deleteLater()
        self.process = None
        if self.config_file and os.path.exists(self.config_file):
            os.remove(self.config_file)
        self.set_running(False)

    def show_crash_report(self) -> None:
        """Show the report which Python wrote to the experiment's log file when the analysis died from a fatal error,
        not otherwise shown in the console."""
        log_files = LOG_FILE_MESSAGE_REGEX.findall(self.log_view.toPlainText())
        if not log_files:
            return
        try:
            log = Path(log_files[-1].strip()).read_text(errors="replace")
        except OSError:
            return
        start = log.rfind(CRASH_REPORT_MARKER)
        if start >= 0:
            self.append_log(f"\n[Crash report, from {log_files[-1].strip()}]\n{log[start:][:MAX_CRASH_REPORT_CHARS]}\n")

    def closeEvent(self, event) -> None:
        if self.is_running():
            reply = QMessageBox.question(
                self, "Analysis running", "An analysis is still running. Stop it and quit?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No, QMessageBox.StandardButton.No,
            )
            if reply != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            process = self.process
            self.stopping = True
            process.terminate()
            if not process.waitForFinished(STOP_TIMEOUT_MS):
                process.kill()
                process.waitForFinished(STOP_TIMEOUT_MS)
            if self.config_file and os.path.exists(self.config_file):
                os.remove(self.config_file)
        event.accept()


def run_config_gui(config: dict | None = None, required_fields: list[str] | None = None) -> None:
    """Show the experiment setup window, from which the analysis can be run, until it is closed."""
    app = QApplication.instance() or QApplication(sys.argv)
    versions = available_versions()
    version = str((config or {}).get("experiment", {}).get("dgame_version", versions[0]))
    defaults, comments = load_defaults(version if version in versions else versions[0])
    model = ConfigModel(defaults, comments)
    if config:
        model.load(config)
    window = ConfigWindow(model, required_fields)
    window.show()
    app.exec()
