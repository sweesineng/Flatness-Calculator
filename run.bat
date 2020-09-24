python.exe -O E:\Python\pyi-mod\Scripts\pyinstaller.exe ^
-y --clean ^
--upx-dir ..\upx-3.96-win64\ ^
--upx-exclude vcruntime140.dll ^
--upx-exclude msvcp140.dll ^
--upx-exclude _uarray.cp38-win_amd64.pyd ^
--add-data "E:\Python\test\flatness-cal-v2.4\icon\app.ico;." ^
--icon="E:\Python\test\flatness-cal-v2.4\icon\app.ico" ^
-F -w ^
-n "Flatness v2.4" ^
flatness-cal-v2.4.py