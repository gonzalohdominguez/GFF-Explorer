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
```
---

## ▶️ Ejecución

Una vez instaladas las dependencias, podés ejecutar el script con:

```bash
python "c:\Users\Ruta_a_a_carpeta\gff_explorer.py"
```
---

## 📁 Formato soportado

El programa trabaja con archivos GFF/GFF3, que contienen 9 columnas estándar:
```bash
seqid | source | type | start | end | score | strand | phase | attributes
```

---

## 🧠 Uso básico

1. Hacer clic en **"Cargar Archivo"**
2. Aplicar filtros desde el panel izquierdo
3. Visualizar y seleccionar registros en la tabla
4. Editar o eliminar filas si es necesario
5. Exportar resultados filtrados
6. Generar estadísticas desde el botón correspondiente

---

## 📊 Estadísticas

El sistema permite generar un gráfico de torta con la distribución de los tipos (`type`), agrupando automáticamente categorías poco frecuentes para mejorar la visualización.

---

## ⚠️ Notas

- Los campos `seqid` y `type` son obligatorios al crear o editar registros.
- Los cambios realizados (edición/eliminación) afectan los datos en memoria.
- Para conservar cambios, es necesario exportar el archivo.

---

## 🛠️ Tecnologías utilizadas

- Python 🐍  
- CustomTkinter 🎨  
- Matplotlib 📊  
- Tkinter (GUI base)

---

## 📄 Licencia

Este proyecto es de uso libre para fines académicos y de investigación.
