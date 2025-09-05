import tkinter as tk
from tkinter import ttk


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
    def __init__(self, root, config):
        self.root = root
        self.config = config
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
        self.result = self.collect_section(self.config, [])
        self.root.destroy()


def initialize_experiment_from_gui(config):
    root = tk.Tk()
    app = ConfigGUI(root, config)
    root.mainloop()
    return app.result
