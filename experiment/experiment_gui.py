import re
import tkinter as tk
from tkinter import filedialog, ttk


class ScrollableFrame(ttk.Frame):
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)

        canvas = tk.Canvas(self)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.scrollable_frame = ttk.Frame(canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")


class ConfigGUI:
    def __init__(self, root, config, required_fields=None):
        self.root = root
        self.config = config
        self.required_fields = required_fields
        self.result = None
        self.entries = {}  # maps path -> (var, type)

        notebook = ttk.Notebook(root)
        notebook.pack(fill="both", expand=True)

        for section_name, section_content in config.items():
            scrollable_tab = ScrollableFrame(notebook)
            notebook.add(scrollable_tab, text=section_name)
            self.build_section(section_content, scrollable_tab.scrollable_frame, path=[section_name])

        submit_btn = ttk.Button(root, text="Run Analysis", command=self.submit)
        submit_btn.pack(pady=10)

    def build_section(self, section, parent, path):
        """Recursively build widgets for config values with path tracking."""
        if isinstance(section, dict):
            for key, value in section.items():
                frame = ttk.LabelFrame(parent, text=key)
                frame.pack(fill="x", padx=5, pady=5, anchor="n")
                self.build_section(value, frame, path + [key])

        elif isinstance(section, list):
            text = tk.Text(parent, height=5, width=40)
            text.insert("1.0", "\n".join(str(v) for v in section))
            text.pack(fill="x", padx=5, pady=5)
            self.entries[".".join(path)] = (text, "list")

        elif isinstance(section, bool):
            var = tk.BooleanVar(value=section)
            cb = ttk.Checkbutton(parent, variable=var, text="enabled")
            cb.pack(anchor="w", padx=5, pady=2)
            self.entries[".".join(path)] = (var, "bool")

        elif isinstance(section, (int, float, str)) or section is None:
            var = tk.StringVar(value="" if section is None else str(section))

            # detect path-like fields
            if any(path[-1].endswith(suffix) for suffix in ["_dir", "_path", "root", "_file"]):
                frame2 = ttk.Frame(parent)
                frame2.pack(fill="x", padx=5, pady=2)

                entry = ttk.Entry(frame2, textvariable=var)
                entry.pack(side="left", fill="x", expand=True)

                def browse(var=var, is_dir=path[-1].endswith("_dir") or path[-1] == "root"):
                    if is_dir:
                        dirname = filedialog.askdirectory()
                        if dirname:
                            var.set(dirname)
                    else:
                        filename = filedialog.askopenfilename()
                        if filename:
                            var.set(filename)

                btn = ttk.Button(frame2, text="Browse", command=browse)
                btn.pack(side="right")

                self.entries[".".join(path)] = (var, "scalar")

            else:
                # regular text entry
                entry = ttk.Entry(parent, textvariable=var)
                entry.pack(fill="x", padx=5, pady=2)
                self.entries[".".join(path)] = (var, "scalar")

    def collect_section(self, section, path):
        """Rebuild nested dict/list/scalar values from entries."""
        if isinstance(section, dict):
            return {k: self.collect_section(v, path + [k]) for k, v in section.items()}

        elif isinstance(section, list):
            widget, wtype = self.entries[".".join(path)]
            text = widget.get("1.0", "end").strip()
            return [v.strip() for v in text.split("\n") if v.strip()]

        elif isinstance(section, bool):
            var, wtype = self.entries[".".join(path)]
            return var.get()

        elif isinstance(section, (int, float, str)) or section is None:
            var, wtype = self.entries[".".join(path)]
            val = var.get()
            if val.isdigit():
                return int(val)
            try:
                return float(val)
            except ValueError:
                return val

        return section

    def submit(self):
        updated = self.collect_section(self.config, [])

        missing = []
        if self.required_fields is not None:
            for field, (var, wtype) in self.entries.items():
                if any(re.search(req_pattern, field) for req_pattern in self.required_fields):
                    if wtype == "scalar":
                        val = var.get().strip()
                        if not val:
                            missing.append(field)

                    elif wtype == "list":
                        text_content = var.get("1.0", "end").strip()
                        if not text_content:
                            missing.append(field)

                    elif wtype == "bool":
                        # booleans are always valid, skip
                        pass

        if missing:
            tk.messagebox.showerror(
                "Missing Required Fields",
                "Please fill in the following required fields:\n" + "\n".join(missing),
            )
            return

        self.result = updated
        self.root.destroy()


def initialize_experiment_from_gui(config: dict,
                                   required_fields: list = None
                                   ) -> dict:
    root = tk.Tk()
    app = ConfigGUI(root, config, required_fields=required_fields)
    root.mainloop()
    return app.result
