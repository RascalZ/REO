# Render only the two M08 document artifacts with installed Microsoft Word.
# Does not open study data, fit models or calculate performance metrics.
$ErrorActionPreference = 'Stop'
$taskRoot = Split-Path -Parent $PSScriptRoot
$qa = Join-Path $taskRoot 'manuscript_v1\qa_m08_v1'
New-Item -ItemType Directory -Force -Path $qa | Out-Null
$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
try {
    $items = @(
        @((Join-Path $taskRoot 'manuscript_v1\Manuscript_BMC_Genomics_M08_v1.docx'), (Join-Path $qa 'main_M08.pdf')),
        @((Join-Path $taskRoot 'manuscript_v1\bmc_additional_files_M08_v1\Additional_file_14.docx'), (Join-Path $qa 'supplement_M08.pdf'))
    )
    foreach ($item in $items) {
        $doc = $null
        try {
            $doc = $word.Documents.Open($item[0], $false, $true)
            $doc.ExportAsFixedFormat($item[1], 17)
            Write-Output $item[1]
        } finally {
            if ($doc) { $doc.Close(0) }
        }
    }
} finally {
    $word.Quit()
    [Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
}
