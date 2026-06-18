import os
import csv
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import customtkinter as ctk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from collections import Counter

# Configuración del tema de la interfaz
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

def parse_comment_region(line):
    """
    Analiza una línea de comentario en busca de información de scaffold/cromosoma (regiones).
    Soporta formatos como:
    ##sequence-region AaegL5_1 1 310827022
    #region AaegL5_1 1 310827022
    # region AaegL5_1 1 310827022
    Retorna un diccionario de registro GFF si coincide, de lo contrario None.
    """
    content = line.lstrip('#').strip()
    tokens = content.split()
    if not tokens:
        return None
        
    idx = -1
    for i, t in enumerate(tokens):
        if t.lower() in ["sequence-region", "region"]:
            idx = i
            break
            
    if idx != -1 and len(tokens) >= idx + 4:
        seqid = tokens[idx + 1]
        start_val = tokens[idx + 2]
        end_val = tokens[idx + 3]
        # Validar que las coordenadas sean numéricas
        if start_val.isdigit() and end_val.isdigit():
            return {
                "seqid": seqid,
                "source": "sequence-region",
                "type": "region",
                "start": start_val,
                "end": end_val,
                "score": ".",
                "strand": ".",
                "phase": ".",
                "attributes": "",
                "gene_id": ""
            }
    return None

def extract_gene_id(attributes):
    """
    Extrae el ID o gene_id del campo de atributos.
    Soporta formato GFF3 (ID=...) y GTF (gene_id "...").
    """
    if not attributes:
        return ""
    
    parts = attributes.split(';')
    for part in parts:
        part = part.strip()
        if not part:
            continue
        
        # Formato GFF3: clave=valor
        if '=' in part:
            k, v = part.split('=', 1)
            k = k.strip()
            v = v.strip()
            if k.upper() in ["ID", "GENE_ID", "GENEID", "GENE"]:
                return v
        
        # Formato GTF: clave "valor"
        else:
            tokens = part.split(None, 1)
            if len(tokens) == 2:
                k, v = tokens
                k = k.strip()
                v = v.strip().strip('"').strip("'")
                if k.upper() in ["GENE_ID", "ID", "GENEID", "GENE"]:
                    return v
    return ""

class GFFExplorerApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        
        self.title("GFF/GFF3 Explorer")
        self.geometry("1300x800")
        
        # Estructuras de datos en memoria y paginación
        self.data = [] # Lista original de diccionarios
        self.filtered_data = [] # Lista de datos filtrados
        self.active_filters = [] # Lista de filtros: {"column": col, "operator": op, "value": val}
        self.current_page = 0
        self.page_size = 5000
        
        self.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        self.setup_ui()
        
    def on_closing(self):
        plt.close('all')
        self.quit()
        self.destroy()
        
    def show_info_dialog(self):
        """Muestra una ventana modal con información institucional y enlace copiable de GitHub."""
        info_window = ctk.CTkToplevel(self)
        info_window.title("Acerca de GFF/GFF3 Explorer")
        info_window.geometry("520x480")
        info_window.resizable(False, False)
        info_window.transient(self)
        info_window.grab_set()
        
        # Encabezado con fondo de acento
        header_frame = ctk.CTkFrame(info_window, fg_color="#1f538d", corner_radius=0)
        header_frame.pack(fill="x")
        
        lbl_icon = ctk.CTkLabel(header_frame, text="🧬", font=ctk.CTkFont(size=36))
        lbl_icon.pack(pady=(18, 0))
        
        lbl_title = ctk.CTkLabel(header_frame, text="GFF/GFF3 Explorer", font=ctk.CTkFont(size=22, weight="bold"), text_color="white")
        lbl_title.pack(pady=(4, 2))
        
        lbl_version = ctk.CTkLabel(header_frame, text="Herramienta de Análisis Genómico · v2.0", font=ctk.CTkFont(size=11), text_color="#b0c9f5")
        lbl_version.pack(pady=(0, 16))
        
        # Cuerpo informativo
        body_frame = ctk.CTkFrame(info_window, fg_color="transparent")
        body_frame.pack(fill="both", expand=True, padx=25, pady=15)
        
        # Institución
        inst_frame = ctk.CTkFrame(body_frame, fg_color="#2b2b2b")
        inst_frame.pack(fill="x", pady=6)
        ctk.CTkLabel(inst_frame, text="🏛  Institución", font=ctk.CTkFont(size=12, weight="bold"), anchor="w").pack(padx=12, pady=(8,2), anchor="w")
        ctk.CTkLabel(inst_frame, text="CENEXA-CREG · UNLP-CONICET", font=ctk.CTkFont(size=12), text_color="#90caf9", anchor="w").pack(padx=12, pady=(0,2), anchor="w")
        ctk.CTkLabel(inst_frame, text="Laboratorio de Neurobiología de Insectos", font=ctk.CTkFont(size=12), anchor="w").pack(padx=12, pady=(0,8), anchor="w")
        
        # Descripción
        desc_frame = ctk.CTkFrame(body_frame, fg_color="#2b2b2b")
        desc_frame.pack(fill="x", pady=6)
        ctk.CTkLabel(desc_frame, text="📄  Descripción", font=ctk.CTkFont(size=12, weight="bold"), anchor="w").pack(padx=12, pady=(8,2), anchor="w")
        ctk.CTkLabel(
            desc_frame,
            text="Exploración, edición, filtrado y análisis de archivos\nde anotación genómica (GFF, GFF3, GTF).",
            font=ctk.CTkFont(size=12), justify="left", anchor="w"
        ).pack(padx=12, pady=(0,8), anchor="w")
        
        # Repositorio GitHub
        git_frame = ctk.CTkFrame(body_frame, fg_color="#2b2b2b")
        git_frame.pack(fill="x", pady=6)
        ctk.CTkLabel(git_frame, text="🐞  Repositorio GitHub", font=ctk.CTkFont(size=12, weight="bold"), anchor="w").pack(padx=12, pady=(8,4), anchor="w")
        
        link_row = ctk.CTkFrame(git_frame, fg_color="transparent")
        link_row.pack(padx=10, pady=(0,10), fill="x")
        
        entry_link = ctk.CTkEntry(link_row, font=ctk.CTkFont(size=11), justify="left")
        entry_link.insert(0, "https://github.com/gonzalohdominguez/GFF-Explorer")
        entry_link.configure(state="readonly")
        entry_link.pack(side="left", fill="x", expand=True, padx=(2, 6))
        
        def copy_link():
            self.clipboard_clear()
            self.clipboard_append("https://github.com/gonzalohdominguez/GFF-Explorer")
            btn_copy.configure(text="✔ Copiado")
            info_window.after(1800, lambda: btn_copy.configure(text="Copiar"))
        
        btn_copy = ctk.CTkButton(link_row, text="Copiar", width=70, command=copy_link)
        btn_copy.pack(side="right", padx=2)
        
        # Pie con derechos
        footer_frame = ctk.CTkFrame(info_window, fg_color="#1a1a1a", corner_radius=0)
        footer_frame.pack(fill="x", side="bottom")
        ctk.CTkLabel(footer_frame, text="Todos los derechos reservados © 2026 · CENEXA-CREG (UNLP-CONICET)", font=ctk.CTkFont(size=9), text_color="gray").pack(pady=8)
        
        # Botón Cerrar
        btn_close = ctk.CTkButton(info_window, text="Cerrar", width=120, command=info_window.destroy)
        btn_close.pack(pady=12, side="bottom")
        
    def setup_ui(self):
        # Configurar la cuadrícula principal (Grid)
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(1, weight=1)
        
        # --- Panel Superior (Barra de Herramientas) ---
        self.top_frame = ctk.CTkFrame(self, height=60)
        self.top_frame.grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=(10,0))
        
        self.btn_load = ctk.CTkButton(self.top_frame, text="Cargar Archivo", command=self.load_file)
        self.btn_load.pack(side="left", padx=10, pady=10)
        
        self.btn_export = ctk.CTkButton(self.top_frame, text="Exportar Datos Filtrados", command=self.export_data)
        self.btn_export.pack(side="left", padx=10, pady=10)
        
        self.btn_stats = ctk.CTkButton(self.top_frame, text="Estadística básica", command=self.show_statistics)
        self.btn_stats.pack(side="left", padx=10, pady=10)
        
        self.btn_info = ctk.CTkButton(self.top_frame, text="Información", fg_color="gray", command=self.show_info_dialog)
        self.btn_info.pack(side="left", padx=10, pady=10)
        
        self.lbl_info = ctk.CTkLabel(self.top_frame, text="Ningún archivo cargado", text_color="gray")
        self.lbl_info.pack(side="right", padx=10, pady=10)
        
        # --- Panel Lateral Izquierdo (Tabview y CRUD) ---
        self.left_frame = ctk.CTkFrame(self, width=280)
        self.left_frame.grid(row=1, column=0, sticky="nswe", padx=10, pady=10)
        self.grid_columnconfigure(0, weight=0) # Barra lateral con ancho fijo
        
        # Tabview para organizar diferentes tipos de filtros (ancho y alto ajustados)
        self.tab_view = ctk.CTkTabview(self.left_frame, width=260, height=300)
        self.tab_view.pack(pady=(5, 5), padx=10, fill="both", expand=False)
        
        self.tab_view.add("Filtros Rápidos")
        self.tab_view.add("Filtro Columna")
        self.tab_view.add("Lista de IDs")
        
        # Pestaña 1: Filtros Rápidos
        tab_rapid = self.tab_view.tab("Filtros Rápidos")
        self.entry_seqid = ctk.CTkEntry(tab_rapid, placeholder_text="seqid...")
        self.entry_seqid.pack(pady=5, padx=5, fill="x")
        
        self.entry_type = ctk.CTkEntry(tab_rapid, placeholder_text="type...")
        self.entry_type.pack(pady=5, padx=5, fill="x")
        
        self.combo_strand = ctk.CTkComboBox(tab_rapid, values=["Todas las hebras", "+", "-", "."])
        self.combo_strand.pack(pady=5, padx=5, fill="x")
        self.combo_strand.set("Todas las hebras")
        
        self.btn_apply_rapid = ctk.CTkButton(tab_rapid, text="Agregar Filtros Rápidos", command=self.apply_rapid_filters)
        self.btn_apply_rapid.pack(pady=10, padx=5, fill="x")
        
        # Pestaña 2: Filtro por Columna
        tab_col = self.tab_view.tab("Filtro Columna")
        self.combo_filter_col = ctk.CTkComboBox(tab_col, values=[
            "seqid", "source", "type", "gene_id", "start", "end", "score", "strand", "phase", "attributes"
        ], command=self.on_filter_column_changed)
        self.combo_filter_col.pack(pady=5, padx=5, fill="x")
        self.combo_filter_col.set("seqid")
        
        self.combo_operator = ctk.CTkComboBox(tab_col, values=[])
        self.combo_operator.pack(pady=5, padx=5, fill="x")
        
        self.entry_filter_val = ctk.CTkEntry(tab_col, placeholder_text="Valor a filtrar...")
        self.entry_filter_val.pack(pady=5, padx=5, fill="x")
        
        self.lbl_value_helper = ctk.CTkLabel(tab_col, text="", text_color="gray", font=ctk.CTkFont(size=10))
        self.lbl_value_helper.pack(pady=2, padx=5, anchor="w")
        
        # Inicializar operadores de forma segura después de que todos los widgets de la pestaña estén creados
        self.on_filter_column_changed("seqid")
        
        self.btn_add_col_filter = ctk.CTkButton(tab_col, text="Agregar Filtro", command=self.add_column_filter)
        self.btn_add_col_filter.pack(pady=10, padx=5, fill="x")
        
        # Pestaña 3: Lista de IDs
        tab_ids = self.tab_view.tab("Lista de IDs")
        self.btn_load_txt = ctk.CTkButton(tab_ids, text="Cargar desde TXT...", command=self.load_ids_from_txt)
        self.btn_load_txt.pack(pady=5, padx=5, fill="x")
        
        lbl_txt_note = ctk.CTkLabel(tab_ids, text="El archivo TXT debe contener un ID por línea.", text_color="gray", font=ctk.CTkFont(size=10), wraplength=200)
        lbl_txt_note.pack(pady=2, padx=5)
        
        self.text_manual_ids = ctk.CTkTextbox(tab_ids, height=100)
        self.text_manual_ids.pack(pady=5, padx=5, fill="both", expand=True)
        self.text_manual_ids.insert("1.0", "Copiar/pegar IDs aquí (separados por comas o saltos de línea)...")
        self.text_manual_ids.bind("<FocusIn>", self.clear_ids_placeholder)
        
        self.btn_apply_ids = ctk.CTkButton(tab_ids, text="Aplicar Lista de IDs", command=self.apply_ids_list)
        self.btn_apply_ids.pack(pady=5, padx=5, fill="x")
        
        # --- Panel de Edición/Gestión (CRUD) Resaltado ---
        self.crud_frame = ctk.CTkFrame(self.left_frame, fg_color="#1f538d", border_width=1, border_color="#1f538d")
        self.crud_frame.pack(pady=(10, 5), padx=10, fill="x")
        
        lbl_crud = ctk.CTkLabel(self.crud_frame, text="📝 GESTIÓN / EDICIÓN", font=ctk.CTkFont(size=12, weight="bold"), text_color="white")
        lbl_crud.pack(pady=(5, 2), padx=10, anchor="w")
        
        buttons_frame = ctk.CTkFrame(self.crud_frame, fg_color="transparent")
        buttons_frame.pack(pady=(2, 5), padx=10, fill="x")
        
        self.btn_add = ctk.CTkButton(buttons_frame, text="Añadir", width=65, fg_color="#2b2b2b", hover_color="#3a3a3a", command=self.add_row)
        self.btn_add.pack(side="left", padx=2, fill="x", expand=True)
        
        self.btn_edit = ctk.CTkButton(buttons_frame, text="Editar", width=65, fg_color="#2b2b2b", hover_color="#3a3a3a", command=self.edit_row)
        self.btn_edit.pack(side="left", padx=2, fill="x", expand=True)
        
        self.btn_delete = ctk.CTkButton(buttons_frame, text="Eliminar", width=65, fg_color="#C62828", hover_color="#B71C1C", command=self.delete_row)
        self.btn_delete.pack(side="left", padx=2, fill="x", expand=True)
        
        # --- Gestión de Filtros Activos ---
        lbl_active_filters = ctk.CTkLabel(self.left_frame, text="Filtros Activos", font=ctk.CTkFont(size=14, weight="bold"))
        lbl_active_filters.pack(pady=(5, 2), padx=10, anchor="w")
        
        self.scrollable_filters_frame = ctk.CTkScrollableFrame(self.left_frame, width=240, height=60)
        self.scrollable_filters_frame.pack(pady=2, padx=10, fill="both", expand=False)
        
        self.btn_clear_filters = ctk.CTkButton(self.left_frame, text="Limpiar Filtros", command=self.clear_filters, fg_color="gray")
        self.btn_clear_filters.pack(pady=5, padx=10, fill="x")
        
        # Firma de desarrollador (en flujo normal para garantizar visibilidad en pantalla)
        lbl_signature = ctk.CTkLabel(self.left_frame, text="Desarrollado por CENEXA-CREG (UNLP-CONICET)", font=ctk.CTkFont(size=9, weight="bold"), text_color="gray")
        lbl_signature.pack(pady=(8, 5), padx=10, fill="x")
        
        # --- Panel Central (Tabla) ---
        self.center_frame = ctk.CTkFrame(self)
        self.center_frame.grid(row=1, column=1, sticky="nswe", padx=(0,10), pady=10)
        
        self.setup_table()
        
    def setup_table(self):
        # Estilo del Treeview
        style = ttk.Style()
        style.theme_use("default")
        style.configure("Treeview", 
                        background="#2b2b2b", 
                        foreground="white", 
                        rowheight=25, 
                        fieldbackground="#2b2b2b")
        style.map("Treeview", background=[("selected", "#1f538d")])
        style.configure("Treeview.Heading", background="#565b5e", foreground="white", font=('Arial', 10, 'bold'))
        
        # 9 columnas estándar + 1 columna auxiliar (gene_id)
        self.columns = ("seqid", "source", "type", "gene_id", "start", "end", "score", "strand", "phase", "attributes")
        
        self.tree = ttk.Treeview(self.center_frame, columns=self.columns, show="headings")
        
        # Vincular Ctrl+c y Ctrl+C para copiar la selección
        self.tree.bind("<Control-c>", self.copy_selection_to_clipboard)
        self.tree.bind("<Control-C>", self.copy_selection_to_clipboard)
        
        vsb = ttk.Scrollbar(self.center_frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(self.center_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        
        self.center_frame.grid_rowconfigure(0, weight=1)
        self.center_frame.grid_columnconfigure(0, weight=1)
        
        widths = {
            "seqid": 100, "source": 110, "type": 90, "gene_id": 120,
            "start": 80, "end": 80, "score": 50, 
            "strand": 50, "phase": 50, "attributes": 350
        }
        for col in self.columns:
            # Asociar el clic en la cabecera con el ordenamiento
            self.tree.heading(col, text=col, command=lambda c=col: self.sort_column(c, False))
            self.tree.column(col, width=widths[col], minwidth=50, stretch=tk.YES if col in ["attributes", "gene_id"] else tk.NO)
            
        # --- Barra de Navegación de Páginas (Paginación) ---
        self.nav_frame = ctk.CTkFrame(self.center_frame, fg_color="transparent")
        self.nav_frame.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(5, 0))
        
        self.btn_prev = ctk.CTkButton(self.nav_frame, text="◀ Anterior", width=100, command=self.prev_page)
        self.btn_prev.pack(side="left", padx=10)
        
        self.lbl_page_info = ctk.CTkLabel(self.nav_frame, text="Página 1 de 1", font=ctk.CTkFont(size=12, weight="bold"))
        self.lbl_page_info.pack(side="left", expand=True)
        
        self.btn_next = ctk.CTkButton(self.nav_frame, text="Siguiente ▶", width=100, command=self.next_page)
        self.btn_next.pack(side="right", padx=10)

    def load_file(self):
        file_path = filedialog.askopenfilename(
            title="Seleccionar Archivo GFF/GFF3",
            filetypes=(("GFF/GFF3 files", "*.gff;*.gff3;*.gtf"), ("All files", "*.*"))
        )
        if not file_path:
            return
            
        self.data.clear()
        parsed_regions_from_comments = []
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    # Capturar y procesar comentarios
                    if line.startswith("#"):
                        reg = parse_comment_region(line)
                        if reg:
                            parsed_regions_from_comments.append(reg)
                        continue
                        
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 9:
                        row_dict = {
                            "seqid": parts[0],
                            "source": parts[1],
                            "type": parts[2],
                            "start": parts[3],
                            "end": parts[4],
                            "score": parts[5],
                            "strand": parts[6],
                            "phase": parts[7],
                            "attributes": parts[8],
                            "gene_id": extract_gene_id(parts[8])
                        }
                        self.data.append(row_dict)
            
            # Evitar duplicados de sequence-region
            existing_regions = {
                (row["seqid"], row["start"], row["end"])
                for row in self.data
                if row["type"].lower() == "region"
            }
            
            for reg in parsed_regions_from_comments:
                key = (reg["seqid"], reg["start"], reg["end"])
                if key not in existing_regions:
                    self.data.append(reg)
                    existing_regions.add(key)
                    
            self.filtered_data = self.data.copy()
            self.current_page = 0 # Reiniciar paginación
            self.active_filters.clear() # Reiniciar filtros activos al cargar un archivo nuevo
            self.refresh_active_filters_display()
            self.refresh_table()
            
        except Exception as e:
            messagebox.showerror("Error de lectura", f"No se pudo procesar el archivo:\n{str(e)}")
            
    def prev_page(self):
        """Retrocede una página en la tabla."""
        if self.current_page > 0:
            self.current_page -= 1
            self.refresh_table()
            
    def next_page(self):
        """Avanza una página en la tabla."""
        total_filtered = len(self.filtered_data)
        max_page = (total_filtered - 1) // self.page_size
        if self.current_page < max_page:
            self.current_page += 1
            self.refresh_table()

    def copy_selection_to_clipboard(self, event):
        """Copia los registros seleccionados en formato tabular (delimitado por tabulaciones) al portapapeles del sistema."""
        selected_items = self.tree.selection()
        if not selected_items:
            return
            
        lines = []
        for item in selected_items:
            values = self.tree.item(item, "values")
            # Unir columnas por tabulaciones para que pegue bien en Excel/Bloc de notas
            lines.append("\t".join(str(val) for val in values))
            
        clipboard_text = "\n".join(lines)
        
        self.clipboard_clear()
        self.clipboard_append(clipboard_text)
        self.update()

    def refresh_table(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
            
        start_idx = self.current_page * self.page_size
        end_idx = start_idx + self.page_size
        
        # Insertar la porción de datos de la página actual
        page_data = self.filtered_data[start_idx:end_idx]
        for idx, row in enumerate(page_data):
            values = tuple(row[col] for col in self.columns)
            self.tree.insert("", "end", iid=str(start_idx + idx), values=values)
            
        # Actualizar información de paginación
        total_filtered = len(self.filtered_data)
        max_page = max(0, (total_filtered - 1) // self.page_size)
        
        self.lbl_page_info.configure(
            text=f"Página {self.current_page + 1} de {max_page + 1} (Registros {start_idx + 1 if total_filtered > 0 else 0:,} - {min(end_idx, total_filtered):,} de {total_filtered:,})"
        )
        
        # Habilitar/Deshabilitar botones de navegación
        if self.current_page == 0:
            self.btn_prev.configure(state="disabled")
        else:
            self.btn_prev.configure(state="normal")
            
        if self.current_page >= max_page:
            self.btn_next.configure(state="disabled")
        else:
            self.btn_next.configure(state="normal")
            
        # Actualizar etiqueta informativa superior
        self.lbl_info.configure(
            text=f"Filtrados: {total_filtered:,} (Total: {len(self.data):,}) | Filtros activos: {len(self.active_filters)}",
            text_color="white"
        )

    def on_filter_column_changed(self, choice):
        # Adaptar operadores de la UI dependiendo de si la columna es numérica o texto
        numeric_cols = ["start", "end", "score", "phase"]
        if choice in numeric_cols:
            self.combo_operator.configure(values=["=", "!=", ">", ">=", "<", "<=", "entre"])
            self.combo_operator.set("=")
            self.lbl_value_helper.configure(text="Para 'entre' use: min,max")
        else:
            self.combo_operator.configure(values=["contiene", "no contiene", "igual a", "distinto de", "comienza con", "termina con"])
            self.combo_operator.set("contiene")
            self.lbl_value_helper.configure(text="")

    def clear_ids_placeholder(self, event):
        current = self.text_manual_ids.get("1.0", "end-1c").strip()
        if "Copiar/pegar IDs aquí" in current:
            self.text_manual_ids.delete("1.0", "end")

    def check_filter_condition(self, row, filt):
        col = filt["column"]
        op = filt["operator"]
        val = filt["value"]
        
        if col not in row:
            return False
            
        row_val = row[col]
        
        # Operador especial para listas
        if op == "en_lista":
            return row_val in val
            
        is_numeric_op = op in ["=", "!=", ">", ">=", "<", "<=", "entre"]
        if is_numeric_op:
            try:
                if op == "entre":
                    parts = [float(x.strip()) for x in val.split(",")]
                    if len(parts) != 2:
                        return False
                    return parts[0] <= float(row_val) <= parts[1]
                else:
                    r_num = float(row_val)
                    f_num = float(val)
                    if op == "=": return r_num == f_num
                    elif op == "!=": return r_num != f_num
                    elif op == ">": return r_num > f_num
                    elif op == ">=": return r_num >= f_num
                    elif op == "<": return r_num < f_num
                    elif op == "<=": return r_num <= f_num
            except ValueError:
                # Si falla el parseo (ej: score es "."), actuar de forma segura
                if op == "!=":
                    return True
                return False
                
        # Comparaciones de cadena
        r_str = str(row_val).strip().lower()
        f_str = str(val).strip().lower()
        
        if op == "contiene":
            return f_str in r_str
        elif op == "no contiene":
            return f_str not in r_str
        elif op == "igual a":
            return r_str == f_str
        elif op == "distinto de":
            return r_str != f_str
        elif op == "comienza con":
            return r_str.startswith(f_str)
        elif op == "termina con":
            return r_str.endswith(f_str)
            
        return False

    def apply_filters(self):
        self.filtered_data = []
        for row in self.data:
            passed = True
            for filt in self.active_filters:
                if not self.check_filter_condition(row, filt):
                    passed = False
                    break
            if passed:
                self.filtered_data.append(row)
                
        self.current_page = 0 # Reiniciar página al cambiar filtros
        self.refresh_table()
        self.refresh_active_filters_display()

    def apply_rapid_filters(self):
        f_seqid = self.entry_seqid.get().strip()
        f_type = self.entry_type.get().strip()
        f_strand = self.combo_strand.get()
        
        any_added = False
        if f_seqid:
            filt = {"column": "seqid", "operator": "contiene", "value": f_seqid}
            if filt not in self.active_filters:
                self.active_filters.append(filt)
                any_added = True
            self.entry_seqid.delete(0, 'end')
            
        if f_type:
            filt = {"column": "type", "operator": "contiene", "value": f_type}
            if filt not in self.active_filters:
                self.active_filters.append(filt)
                any_added = True
            self.entry_type.delete(0, 'end')
            
        if f_strand != "Todas las hebras":
            filt = {"column": "strand", "operator": "igual a", "value": f_strand}
            if filt not in self.active_filters:
                self.active_filters.append(filt)
                any_added = True
            self.combo_strand.set("Todas las hebras")
            
        if any_added:
            self.apply_filters()

    def add_column_filter(self):
        col = self.combo_filter_col.get()
        op = self.combo_operator.get()
        val = self.entry_filter_val.get().strip()
        
        if not val:
            messagebox.showwarning("Atención", "Por favor ingresa un valor para filtrar.")
            return
            
        filt = {"column": col, "operator": op, "value": val}
        if filt not in self.active_filters:
            self.active_filters.append(filt)
            self.apply_filters()
            
        self.entry_filter_val.delete(0, 'end')

    def load_ids_from_txt(self):
        messagebox.showinfo("Aclaración", "El archivo TXT debe contener un ID por línea.")
        file_path = filedialog.askopenfilename(
            title="Seleccionar Archivo de IDs (TXT)",
            filetypes=(("Text files", "*.txt"), ("All files", "*.*"))
        )
        if not file_path:
            return
            
        try:
            ids = set()
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    cleaned = line.strip()
                    if cleaned:
                        ids.add(cleaned)
                        
            if ids:
                filt = {
                    "column": "gene_id",
                    "operator": "en_lista",
                    "value": ids,
                    "filename": os.path.basename(file_path)
                }
                # Reemplazar cualquier filtro anterior por lista de IDs para evitar solapamientos
                self.active_filters = [f for f in self.active_filters if f["operator"] != "en_lista"]
                self.active_filters.append(filt)
                self.apply_filters()
                
        except Exception as e:
            messagebox.showerror("Error", f"No se pudo leer el archivo de IDs:\n{str(e)}")

    def apply_ids_list(self):
        input_text = self.text_manual_ids.get("1.0", "end-1c").strip()
        if not input_text or "Copiar/pegar IDs aquí" in input_text:
            messagebox.showwarning("Atención", "Por favor ingresa o pega algunos IDs.")
            return
            
        normalized = input_text.replace(',', '\n')
        lines = normalized.split('\n')
        ids = set()
        for line in lines:
            cleaned = line.strip()
            if cleaned:
                ids.add(cleaned)
                
        if ids:
            filt = {
                "column": "gene_id",
                "operator": "en_lista",
                "value": ids,
                "filename": "Lista manual"
            }
            # Reemplazar cualquier filtro anterior por lista de IDs
            self.active_filters = [f for f in self.active_filters if f["operator"] != "en_lista"]
            self.active_filters.append(filt)
            self.apply_filters()
            
            # Limpiar text box y colocar placeholder
            self.text_manual_ids.delete("1.0", "end")
            self.text_manual_ids.insert("1.0", "Copiar/pegar IDs aquí (separados por comas o saltos de línea)...")

    def remove_single_filter(self, filter_idx):
        if 0 <= filter_idx < len(self.active_filters):
            self.active_filters.pop(filter_idx)
            self.apply_filters()

    def clear_filters(self):
        self.active_filters.clear()
        self.current_page = 0 # Reiniciar paginación
        
        # Limpiar entradas visuales
        self.entry_seqid.delete(0, 'end')
        self.entry_type.delete(0, 'end')
        self.combo_strand.set("Todas las hebras")
        self.entry_filter_val.delete(0, 'end')
        
        self.text_manual_ids.delete("1.0", "end")
        self.text_manual_ids.insert("1.0", "Copiar/pegar IDs aquí (separados por comas o saltos de línea)...")
        
        self.apply_filters()

    def refresh_active_filters_display(self):
        for widget in self.scrollable_filters_frame.winfo_children():
            widget.destroy()
            
        for idx, filt in enumerate(self.active_filters):
            item_frame = ctk.CTkFrame(self.scrollable_filters_frame, fg_color="#333333")
            item_frame.pack(fill="x", pady=2, padx=2)
            
            col = filt["column"]
            op = filt["operator"]
            val = filt["value"]
            
            if op == "en_lista":
                fname = filt.get("filename", "Lista")
                desc = f"gene_id ∈ {fname} ({len(val)} IDs)"
            elif op in ["=", "!=", ">", ">=", "<", "<=", "entre"]:
                desc = f"{col} {op} {val}"
            else:
                desc = f"{col} {op} '{val}'"
                
            lbl = ctk.CTkLabel(item_frame, text=desc, font=ctk.CTkFont(size=11), anchor="w")
            lbl.pack(side="left", padx=5, fill="x", expand=True)
            
            btn_del = ctk.CTkButton(
                item_frame, 
                text="x", 
                width=18, 
                height=18, 
                fg_color="#C62828", 
                hover_color="#B71C1C",
                command=lambda i=idx: self.remove_single_filter(i)
            )
            btn_del.pack(side="right", padx=5)

    def get_selected_filtered_index(self):
        selected = self.tree.selection()
        if not selected:
            messagebox.showwarning("Atención", "Por favor, selecciona una fila primero.")
            return None
        return int(selected[0])

    def show_crud_dialog(self):
        """Abre una ventana flotante con los controles de edición CRUD."""
        dlg = ctk.CTkToplevel(self)
        dlg.title("Gestión / Edición")
        dlg.geometry("320x200")
        dlg.resizable(False, False)
        dlg.transient(self)
        dlg.grab_set()

        header = ctk.CTkFrame(dlg, fg_color="#7b3f00", corner_radius=0)
        header.pack(fill="x")
        ctk.CTkLabel(header, text="✏️  Gestión / Edición de Registros",
                     font=ctk.CTkFont(size=13, weight="bold"), text_color="white").pack(pady=12, padx=14)

        body = ctk.CTkFrame(dlg, fg_color="transparent")
        body.pack(fill="both", expand=True, padx=20, pady=16)

        def run_and_close(fn):
            dlg.destroy()
            fn()

        self.btn_add = ctk.CTkButton(body, text="➕  Añadir Fila",
                                     fg_color="#2b5f2b", hover_color="#3a7a3a",
                                     font=ctk.CTkFont(size=12),
                                     command=lambda: run_and_close(self.add_row))
        self.btn_add.pack(fill="x", pady=4)

        self.btn_edit = ctk.CTkButton(body, text="🖊  Editar Fila Seleccionada",
                                      fg_color="#2b2b5f", hover_color="#3a3a7a",
                                      font=ctk.CTkFont(size=12),
                                      command=lambda: run_and_close(self.edit_row))
        self.btn_edit.pack(fill="x", pady=4)

        self.btn_delete = ctk.CTkButton(body, text="🗑  Eliminar Fila Seleccionada",
                                        fg_color="#C62828", hover_color="#B71C1C",
                                        font=ctk.CTkFont(size=12),
                                        command=lambda: run_and_close(self.delete_row))
        self.btn_delete.pack(fill="x", pady=4)

    def add_row(self):
        self.open_edit_dialog(mode="add")

    def edit_row(self):
        idx = self.get_selected_filtered_index()
        if idx is not None:
            self.open_edit_dialog(mode="edit", index=idx)

    def delete_row(self):
        idx = self.get_selected_filtered_index()
        if idx is not None:
            confirm = messagebox.askyesno("Confirmar Eliminación", "¿Estás seguro de que deseas eliminar la fila seleccionada de forma permanente?")
            if confirm:
                row_to_delete = self.filtered_data[idx]
                if row_to_delete in self.data:
                    self.data.remove(row_to_delete)
                self.apply_filters()

    def open_edit_dialog(self, mode="add", index=None):
        dialog = ctk.CTkToplevel(self)
        dialog.title("Añadir Nueva Fila" if mode == "add" else "Editar Fila")
        dialog.geometry("400x550")
        dialog.transient(self)
        dialog.grab_set()
        
        entries = {}
        # Únicamente las 9 columnas originales se editan
        gff_cols = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
        for i, col in enumerate(gff_cols):
            lbl = ctk.CTkLabel(dialog, text=col)
            lbl.grid(row=i, column=0, padx=15, pady=10, sticky="e")
            entry = ctk.CTkEntry(dialog, width=250)
            entry.grid(row=i, column=1, padx=10, pady=10, sticky="w")
            entries[col] = entry
            
        if mode == "edit":
            row_data = self.filtered_data[index]
            for col in gff_cols:
                entries[col].insert(0, row_data[col])
                
        def save_data():
            new_row = {}
            for col in gff_cols:
                new_row[col] = entries[col].get()
                
            # Extraer y asignar el gene_id automáticamente
            new_row["gene_id"] = extract_gene_id(new_row["attributes"])
            
            if not new_row["seqid"] or not new_row["type"]:
                messagebox.showerror("Error de Validación", "Los campos 'seqid' y 'type' son obligatorios.", parent=dialog)
                return
                
            if mode == "add":
                self.data.append(new_row)
            elif mode == "edit":
                orig_row = self.filtered_data[index]
                if orig_row in self.data:
                    orig_idx = self.data.index(orig_row)
                    self.data[orig_idx] = new_row
                    
            self.apply_filters()
            dialog.destroy()
            
        btn_save = ctk.CTkButton(dialog, text="Guardar Cambios", command=save_data)
        btn_save.grid(row=len(gff_cols), column=0, columnspan=2, pady=20)

    def export_data(self):
        if not self.filtered_data:
            messagebox.showwarning("Atención", "No hay datos filtrados para exportar.")
            return
            
        file_path = filedialog.asksaveasfilename(
            title="Exportar Datos Filtrados",
            defaultextension=".csv",
            filetypes=(("CSV (Delimitado por comas)", "*.csv"), ("GTF files", "*.gtf"), ("GFF files", "*.gff"), ("All files", "*.*"))
        )
        
        if not file_path:
            return
            
        # Exportar estrictamente las 9 columnas originales
        export_cols = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
        
        try:
            if file_path.endswith('.csv'):
                with open(file_path, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.DictWriter(f, fieldnames=export_cols)
                    writer.writeheader()
                    for row in self.filtered_data:
                        # Crear corte de diccionario excluyendo gene_id
                        row_to_write = {col: row[col] for col in export_cols}
                        writer.writerow(row_to_write)
            else:
                with open(file_path, 'w', encoding='utf-8') as f:
                    for row in self.filtered_data:
                        attrs = row['attributes']
                        if file_path.endswith('.gtf'):
                            new_attrs = []
                            for pair in attrs.split(';'):
                                if '=' in pair:
                                    k, v = pair.split('=', 1)
                                    new_attrs.append(f'{k} "{v}"')
                                else:
                                    if pair: new_attrs.append(pair)
                            attrs = "; ".join(new_attrs) + ";" if new_attrs else attrs
                            
                        line = f"{row['seqid']}\t{row['source']}\t{row['type']}\t{row['start']}\t{row['end']}\t{row['score']}\t{row['strand']}\t{row['phase']}\t{attrs}\n"
                        f.write(line)
                        
            messagebox.showinfo("Exportación Exitosa", f"Se exportaron {len(self.filtered_data)} registros a:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Error de Exportación", str(e))

    def sort_column(self, col, reverse):
        """Ordena los datos filtrados en memoria de manera bidireccional al hacer clic en la cabecera."""
        is_numeric = col in ["start", "end", "score", "phase"]
        
        def get_sort_key(row):
            val = row[col]
            if is_numeric:
                try:
                    return (0, float(val))
                except ValueError:
                    return (1, str(val))
            return (0, str(val).lower())
            
        self.filtered_data.sort(key=get_sort_key, reverse=reverse)
        
        # Quitar indicadores de ordenación de otras columnas
        for c in self.columns:
            self.tree.heading(c, text=c)
            
        # Añadir indicador visual (flecha) a la columna actualmente ordenada
        arrow = " ▼" if reverse else " ▲"
        self.tree.heading(col, text=col + arrow, command=lambda: self.sort_column(col, not reverse))
        
        self.refresh_table()

    def show_statistics(self):
        """Muestra una ventana modal con estadísticas básicas y gráficos de distribución y frecuencia."""
        if not self.filtered_data:
            messagebox.showwarning("Atención", "No hay datos para generar estadísticas.")
            return
            
        total_records = len(self.filtered_data)
        
        # Identificar scaffolds y genes únicos
        scaffolds = set(row['seqid'] for row in self.filtered_data)
        num_scaffolds = len(scaffolds)
        
        genes = set(row['gene_id'] for row in self.filtered_data if row['gene_id'])
        num_genes = len(genes)
        
        # Frecuencia de tipos
        type_counts = Counter(row['type'] for row in self.filtered_data)
        num_exons = type_counts.get("exon", 0)
        num_mrna = type_counts.get("mRNA", 0) or type_counts.get("transcript", 0)
        num_cds = type_counts.get("CDS", 0)
        num_regions = type_counts.get("region", 0)
        
        # Hebras
        strand_counts = Counter(row['strand'] for row in self.filtered_data)
        plus_strand = strand_counts.get("+", 0)
        minus_strand = strand_counts.get("-", 0)
        
        stats_window = ctk.CTkToplevel(self)
        stats_window.title("Estadística Básica del GFF/GTF")
        stats_window.geometry("1000x600")
        stats_window.transient(self)
        stats_window.grab_set()
        
        # Rejilla principal
        stats_window.grid_columnconfigure(0, weight=1)
        stats_window.grid_columnconfigure(1, weight=2)
        stats_window.grid_rowconfigure(0, weight=1)
        
        # Resumen métrico (Izquierda)
        summary_frame = ctk.CTkFrame(stats_window, width=280)
        summary_frame.grid(row=0, column=0, sticky="nswe", padx=15, pady=15)
        
        lbl_title = ctk.CTkLabel(summary_frame, text="Resumen de Métricas", font=ctk.CTkFont(size=18, weight="bold"))
        lbl_title.pack(pady=15, padx=15, anchor="w")
        
        metrics = [
            ("Total Registros:", f"{total_records:,}"),
            ("Scaffolds (seqid):", f"{num_scaffolds:,}"),
            ("Genes Únicos (attributes ID):", f"{num_genes:,}"),
            ("ARNm / Transcritos:", f"{num_mrna:,}"),
            ("Exones:", f"{num_exons:,}"),
            ("CDS:", f"{num_cds:,}"),
            ("Regiones genomic:", f"{num_regions:,}"),
            ("Hebra (+):", f"{plus_strand:,}"),
            ("Hebra (-):", f"{minus_strand:,}")
        ]
        
        for label, val in metrics:
            metric_frame = ctk.CTkFrame(summary_frame, fg_color="transparent")
            metric_frame.pack(fill="x", padx=15, pady=5)
            
            lbl_name = ctk.CTkLabel(metric_frame, text=label, font=ctk.CTkFont(size=12, weight="bold"), anchor="w")
            lbl_name.pack(side="left")
            
            lbl_val = ctk.CTkLabel(metric_frame, text=val, font=ctk.CTkFont(size=12), anchor="e")
            lbl_val.pack(side="right")
            
        # Panel Derecho (Gráficos)
        plots_frame = ctk.CTkFrame(stats_window)
        plots_frame.grid(row=0, column=1, sticky="nswe", padx=15, pady=15)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
        fig.patch.set_facecolor('#2b2b2b')
        
        # Gráfico 1: Torta (Pie Chart) de tipos principales (Top 5 + Otros)
        top_types = type_counts.most_common(5)
        top_labels = [k for k, _ in top_types]
        top_sizes = [v for _, v in top_types]
        other_sum = sum(type_counts.values()) - sum(top_sizes)
        
        if other_sum > 0:
            top_labels.append("Otros")
            top_sizes.append(other_sum)
            
        ax1.set_facecolor('#2b2b2b')
        ax1.pie(top_sizes, labels=top_labels, autopct='%1.1f%%', startangle=140, textprops={'color':"w", 'size': 8})
        ax1.set_title("Distribución de Elementos", color="w", fontsize=11, weight="bold")
        ax1.axis('equal')
        
        # Gráfico 2: Barra Horizontal del Top 5 Scaffolds con más anotaciones
        scaffold_counts = Counter(row['seqid'] for row in self.filtered_data)
        top_scaffolds = scaffold_counts.most_common(5)
        scaf_labels = [k for k, _ in top_scaffolds]
        scaf_sizes = [v for _, v in top_scaffolds]
        
        ax2.set_facecolor('#2b2b2b')
        ax2.tick_params(colors='w', labelsize=8)
        ax2.spines['bottom'].set_color('w')
        ax2.spines['left'].set_color('w')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        
        bars = ax2.bar(scaf_labels, scaf_sizes, color="#1f538d")
        ax2.set_title("Top 5 Scaffolds con Más Elementos", color="w", fontsize=11, weight="bold")
        ax2.set_ylabel("Frecuencia", color="w", fontsize=9)
        plt.setp(ax2.get_xticklabels(), rotation=30, ha="right")
        
        for bar in bars:
            height = bar.get_height()
            ax2.annotate(f'{height:,}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', color='w', fontsize=7)
        
        fig.tight_layout()
        
        canvas = FigureCanvasTkAgg(fig, master=plots_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)

if __name__ == "__main__":
    app = GFFExplorerApp()
    app.mainloop()
