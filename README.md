# Flatness-Calculator
Calculate flatness of surface from point clouds in python gui

## How-to compile to exe
Create a new venv name pyi-env\
  ```python.exe -m venv E:\Python\pyi-env```\
  
Install all dependencies\
  ```pip install numpy pandas scipy matplotlib pyinstaller```\
Note: pyinstaller had some problem with latest matplotlib due to the mpl-data.\
Just change "\pyi-env\Lib\site-packages\PyInstaller\hooks\hook-matplotlib.py" to:
>datas = [
    (mpl_data_dir, "matplotlib/mpl-data"),
]

Download latest UPX and extract it\
Edit the run.bat to change:
- upx-dir
- name
- filename

Excute the run.bat\
```.\run.bat```
