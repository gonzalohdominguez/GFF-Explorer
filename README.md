# 🧬 GFF/GFF3 Explorer - Bioinformatics GUI

Aplicación de escritorio desarrollada en Python para explorar, filtrar, editar y visualizar archivos **GFF/GFF3** de forma interactiva mediante una interfaz gráfica basada en `CustomTkinter`.

---

## 🚀 Características

- 📂 Carga de archivos **GFF/GFF3**
- 🔎 Filtros interactivos por:
  - `seqid`
  - `type`
  - `strand`
- 📊 Visualización de estadísticas (gráfico de torta)
- ✏️ Edición de filas (CRUD completo):
  - Añadir
  - Editar
  - Eliminar
- 💾 Exportación de datos filtrados a:
  - CSV
  - GFF
  - GTF

---

## 📦 Requisitos

Antes de ejecutar el programa, asegurate de tener instaladas las dependencias necesarias.

Abrí tu terminal de **PowerShell** y ejecutá:

```bash
pip install customtkinter matplotlib

---

## ▶️ Ejecución

Una vez instaladas las dependencias, podés ejecutar el script con:

```bash
python "c:\Users\Ruta_a_a_carpeta\gff_explorer.py"

---

## 📁 Formato soportado

El programa trabaja con archivos GFF/GFF3, que contienen 9 columnas estándar:
```bash
seqid | source | type | start | end | score | strand | phase | attributes

---

## 📊 Estadísticas
El sistema permite generar un gráfico de torta con la distribución de los tipos (type), agrupando automáticamente categorías poco frecuentes para mejorar la visualización.

---

📄 Licencia

Este proyecto es de uso libre para fines académicos y de investigación.
