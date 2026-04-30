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

class GFFExplorerApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        
        self.title("GFF/GFF3 Explorer - Bioinformática GUI")
        self.geometry("1250x750")
        
        # Estructuras de datos en memoria
        self.data = [] # Lista original de diccionarios (1 por fila GFF)
        self.filtered_data = [] # Lista que contiene los datos filtrados actuales
        
        self.setup_ui()
        
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
        
        self.btn_stats = ctk.CTkButton(self.top_frame, text="Estadísticas de Tipos", command=self.show_statistics)
        self.btn_stats.pack(side="left", padx=10, pady=10)
        
        self.lbl_info = ctk.CTkLabel(self.top_frame, text="Ningún archivo cargado", text_color="gray")
        self.lbl_info.pack(side="right", padx=10, pady=10)
        
        # --- Panel Lateral Izquierdo (Filtros y CRUD) ---
        self.left_frame = ctk.CTkFrame(self, width=250)
        self.left_frame.grid(row=1, column=0, sticky="nswe", padx=10, pady=10)
        self.left_frame.grid_propagate(False) # Mantener el ancho
        
        # Sección Filtros
        lbl_filters = ctk.CTkLabel(self.left_frame, text="Filtros Interactivos", font=ctk.CTkFont(size=16, weight="bold"))
        lbl_filters.pack(pady=(15, 10), padx=10, anchor="w")
        
        self.entry_seqid = ctk.CTkEntry(self.left_frame, placeholder_text="Filtrar seqid...")
        self.entry_seqid.pack(pady=5, padx=10, fill="x")
        
        self.entry_type = ctk.CTkEntry(self.left_frame, placeholder_text="Filtrar type...")
        self.entry_type.pack(pady=5, padx=10, fill="x")
        
        self.combo_strand = ctk.CTkComboBox(self.left_frame, values=["Todas las hebras", "+", "-", "."])
        self.combo_strand.pack(pady=5, padx=10, fill="x")
        self.combo_strand.set("Todas las hebras")
        
        self.btn_apply_filters = ctk.CTkButton(self.left_frame, text="Aplicar Filtros", command=self.apply_filters)
        self.btn_apply_filters.pack(pady=10, padx=10, fill="x")
        
        self.btn_clear_filters = ctk.CTkButton(self.left_frame, text="Limpiar Filtros", command=self.clear_filters, fg_color="gray")
        self.btn_clear_filters.pack(pady=5, padx=10, fill="x")
        
        # Separador
        ctk.CTkLabel(self.left_frame, text="─"*25, text_color="gray").pack(pady=10)
        
        # Sección CRUD
        lbl_crud = ctk.CTkLabel(self.left_frame, text="Gestión de Filas", font=ctk.CTkFont(size=16, weight="bold"))
        lbl_crud.pack(pady=(5, 10), padx=10, anchor="w")
        
        self.btn_add = ctk.CTkButton(self.left_frame, text="Añadir Nueva Fila", command=self.add_row)
        self.btn_add.pack(pady=5, padx=10, fill="x")
        
        self.btn_edit = ctk.CTkButton(self.left_frame, text="Editar Seleccionada", command=self.edit_row)
        self.btn_edit.pack(pady=5, padx=10, fill="x")
        
        self.btn_delete = ctk.CTkButton(self.left_frame, text="Eliminar Seleccionada", command=self.delete_row, fg_color="#C62828", hover_color="#B71C1C")
        self.btn_delete.pack(pady=5, padx=10, fill="x")
        
        # --- Panel Central (Tabla) ---
        self.center_frame = ctk.CTkFrame(self)
        self.center_frame.grid(row=1, column=1, sticky="nswe", padx=(0,10), pady=10)
        
        self.setup_table()
        
    def setup_table(self):
        # Configurar el estilo del Treeview para que combine con el tema oscuro
        style = ttk.Style()
        style.theme_use("default")
        style.configure("Treeview", 
                        background="#2b2b2b", 
                        foreground="white", 
                        rowheight=25, 
                        fieldbackground="#2b2b2b")
        style.map("Treeview", background=[("selected", "#1f538d")])
        style.configure("Treeview.Heading", background="#565b5e", foreground="white", font=('Arial', 10, 'bold'))
        
        # 9 columnas del formato GFF
        self.columns = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
        
        # Crear Treeview
        self.tree = ttk.Treeview(self.center_frame, columns=self.columns, show="headings")
        
        # Barras de desplazamiento
        vsb = ttk.Scrollbar(self.center_frame, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(self.center_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        
        # Ubicar elementos en el grid del frame central
        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        
        self.center_frame.grid_rowconfigure(0, weight=1)
        self.center_frame.grid_columnconfigure(0, weight=1)
        
        # Definir cabeceras y anchos
        widths = {
            "seqid": 100, "source": 100, "type": 100, 
            "start": 80, "end": 80, "score": 50, 
            "strand": 50, "phase": 50, "attributes": 400
        }
        for col in self.columns:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=widths[col], minwidth=50, stretch=tk.YES if col == "attributes" else tk.NO)

    def load_file(self):
        file_path = filedialog.askopenfilename(
            title="Seleccionar Archivo GFF/GFF3",
            filetypes=(("GFF files", "*.gff;*.gff3"), ("All files", "*.*"))
        )
        if not file_path:
            return
            
        self.data.clear()
        
        try:
            # Lectura y manejo de excepciones
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    # Ignorar comentarios del encabezado
                    if line.startswith("#"):
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
                            "attributes": parts[8]
                        }
                        self.data.append(row_dict)
                        
            # Inicialmente, los datos filtrados son todos los datos
            self.filtered_data = self.data.copy()
            self.lbl_info.configure(text=f"Cargados {len(self.data)} registros", text_color="white")
            self.refresh_table()
            
        except Exception as e:
            messagebox.showerror("Error de lectura", f"No se pudo procesar el archivo:\n{str(e)}")
            
    def refresh_table(self):
        # Limpiar filas existentes
        for item in self.tree.get_children():
            self.tree.delete(item)
            
        # Insertar nuevos datos usando el índice como IID para facilitar operaciones CRUD
        for idx, row in enumerate(self.filtered_data):
            values = tuple(row[col] for col in self.columns)
            self.tree.insert("", "end", iid=str(idx), values=values)

    def apply_filters(self):
        f_seqid = self.entry_seqid.get().strip().lower()
        f_type = self.entry_type.get().strip().lower()
        f_strand = self.combo_strand.get()
        
        self.filtered_data = []
        for row in self.data:
            match_seqid = f_seqid in row['seqid'].lower() if f_seqid else True
            match_type = f_type in row['type'].lower() if f_type else True
            match_strand = (row['strand'] == f_strand) if f_strand != "Todas las hebras" else True
            
            if match_seqid and match_type and match_strand:
                self.filtered_data.append(row)
                
        self.refresh_table()
        self.lbl_info.configure(text=f"Mostrando {len(self.filtered_data)} de {len(self.data)} registros")

    def clear_filters(self):
        self.entry_seqid.delete(0, 'end')
        self.entry_type.delete(0, 'end')
        self.combo_strand.set("Todas las hebras")
        self.filtered_data = self.data.copy()
        self.refresh_table()
        self.lbl_info.configure(text=f"Mostrando {len(self.filtered_data)} de {len(self.data)} registros")

    def get_selected_filtered_index(self):
        selected = self.tree.selection()
        if not selected:
            messagebox.showwarning("Atención", "Por favor, selecciona una fila primero.")
            return None
        return int(selected[0])

    def add_row(self):
        self.open_edit_dialog(mode="add")

    def edit_row(self):
        idx = self.get_selected_filtered_index()
        if idx is not None:
            self.open_edit_dialog(mode="edit", index=idx)

    def delete_row(self):
        idx = self.get_selected_filtered_index()
        if idx is not None:
            # Confirmación explícita
            confirm = messagebox.askyesno("Confirmar Eliminación", "¿Estás seguro de que deseas eliminar la fila seleccionada de forma permanente?")
            if confirm:
                row_to_delete = self.filtered_data[idx]
                
                # Remover de los datos originales si existe
                if row_to_delete in self.data:
                    self.data.remove(row_to_delete)
                
                # Refrescar la vista actual (reaplicar filtros para actualizar `self.filtered_data` y el Treeview)
                self.apply_filters()

    def open_edit_dialog(self, mode="add", index=None):
        dialog = ctk.CTkToplevel(self)
        dialog.title("Añadir Nueva Fila" if mode == "add" else "Editar Fila")
        dialog.geometry("400x550")
        dialog.transient(self) # Hacer que se mantenga por encima de la ventana principal
        dialog.grab_set() # Comportamiento modal
        
        entries = {}
        for i, col in enumerate(self.columns):
            lbl = ctk.CTkLabel(dialog, text=col)
            lbl.grid(row=i, column=0, padx=15, pady=10, sticky="e")
            entry = ctk.CTkEntry(dialog, width=250)
            entry.grid(row=i, column=1, padx=10, pady=10, sticky="w")
            entries[col] = entry
            
        if mode == "edit":
            row_data = self.filtered_data[index]
            for col in self.columns:
                entries[col].insert(0, row_data[col])
                
        def save_data():
            new_row = {}
            for col in self.columns:
                new_row[col] = entries[col].get()
                
            # Validaciones básicas
            if not new_row["seqid"] or not new_row["type"]:
                messagebox.showerror("Error de Validación", "Los campos 'seqid' y 'type' son obligatorios.", parent=dialog)
                return
                
            if mode == "add":
                self.data.append(new_row)
            elif mode == "edit":
                # Buscar en self.data la referencia y modificarla
                orig_row = self.filtered_data[index]
                if orig_row in self.data:
                    orig_idx = self.data.index(orig_row)
                    self.data[orig_idx] = new_row
                    
            self.apply_filters()
            dialog.destroy()
            
        btn_save = ctk.CTkButton(dialog, text="Guardar Cambios", command=save_data)
        btn_save.grid(row=len(self.columns), column=0, columnspan=2, pady=20)

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
            
        try:
            if file_path.endswith('.csv'):
                with open(file_path, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.DictWriter(f, fieldnames=self.columns)
                    writer.writeheader()
                    writer.writerows(self.filtered_data)
            else: # GTF o GFF
                with open(file_path, 'w', encoding='utf-8') as f:
                    for row in self.filtered_data:
                        attrs = row['attributes']
                        # Si es .gtf, intentar conversión básica
                        if file_path.endswith('.gtf'):
                            new_attrs = []
                            # GFF usa ; y =. GTF usa ; y espacio con comillas.
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

    def show_statistics(self):
        if not self.filtered_data:
            messagebox.showwarning("Atención", "No hay datos para generar estadísticas.")
            return
            
        # Calcular frecuencias
        type_counts = Counter(row['type'] for row in self.filtered_data)
        
        if not type_counts:
            return
            
        # Para que el gráfico no se vuelva ilegible con muchos tipos, agrupar los menores al 2% en 'Otros' (Opcional, pero recomendado en Bioinformática)
        total = sum(type_counts.values())
        labels = []
        sizes = []
        otros = 0
        
        for k, v in type_counts.items():
            if (v / total) < 0.015 and len(type_counts) > 10: # Agrupar si es menor al 1.5% y hay muchos
                otros += v
            else:
                labels.append(k)
                sizes.append(v)
                
        if otros > 0:
            labels.append('Otros')
            sizes.append(otros)
        
        # Ventana para el gráfico
        stats_window = ctk.CTkToplevel(self)
        stats_window.title("Frecuencia de Tipos (Type)")
        stats_window.geometry("700x550")
        
        # Figura matplotlib
        fig, ax = plt.subplots(figsize=(7, 5))
        
        # Colores consistentes y modo oscuro en plot
        fig.patch.set_facecolor('#2b2b2b')
        ax.set_facecolor('#2b2b2b')
        
        wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, textprops={'color':"w"})
        
        # Ajuste de layout para que no se corten las etiquetas
        ax.axis('equal')  
        plt.setp(autotexts, size=8, weight="bold")
        
        # Embeber matplotlib en tkinter
        canvas = FigureCanvasTkAgg(fig, master=stats_window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=20, pady=20)

if __name__ == "__main__":
    app = GFFExplorerApp()
    app.mainloop()
