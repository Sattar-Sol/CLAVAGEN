@echo off
REM Creates an "CLAVAGEN" icon on the Desktop that starts the program.
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -Command ^
 "$s=(New-Object -ComObject WScript.Shell).CreateShortcut([Environment]::GetFolderPath('Desktop')+'\CLAVAGEN.lnk');" ^
 "$s.TargetPath='%~dp0Start_Windows.bat'; $s.WorkingDirectory='%~dp0'; $s.IconLocation='%~dp0tools\AppIcon.ico';" ^
 "$s.WindowStyle=7; $s.Description='CLAVAGEN clinical variant calling'; $s.Save()"
echo Desktop shortcut "CLAVAGEN" created.
pause
