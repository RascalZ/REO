# Local author reading-copy export; never opens the legacy manuscript or study data.
$ErrorActionPreference = 'Stop'
# VBScript dispatch avoids the observed PowerShell COM export stall on this host.
& cscript.exe //nologo (Join-Path $PSScriptRoot 'render_author_review_v1.vbs')
if ($LASTEXITCODE -ne 0) { throw 'Word export failed' }
