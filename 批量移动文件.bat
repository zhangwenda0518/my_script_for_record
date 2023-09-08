@echo off
for /f "delims=" %%a in ('dir /a-d /b /s ') do (move "%%~a" ./)