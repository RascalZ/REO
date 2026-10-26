' Read-only export of the two revised manuscript reading copies.
Option Explicit
Dim fso, taskRoot, app, doc, relativeBase, fileBase, freshPdf, finalPdf
Set fso = CreateObject("Scripting.FileSystemObject")
taskRoot = fso.GetParentFolderName(fso.GetParentFolderName(WScript.ScriptFullName))
Set app = CreateObject("Word.Application")
app.Visible = False
app.DisplayAlerts = 0
For Each relativeBase In Array("manuscript_v1\author_review_v3\Manuscript_BMC_Genomics_author_review_v3", "manuscript_v1\author_review_v3_zh\Manuscript_BMC_Genomics_author_review_v3_zh")
    fileBase = fso.BuildPath(taskRoot, relativeBase)
    freshPdf = fileBase & ".vb_rendered.pdf"
    finalPdf = fileBase & ".pdf"
    If fso.FileExists(freshPdf) Or fso.FileExists(finalPdf) Then Err.Raise 513, , "Preserve existing PDF copies"
    WScript.Echo "Opening " & relativeBase
    Set doc = app.Documents.Open(fileBase & ".docx", False, True)
    doc.Repaginate
    doc.Fields.Update
    WScript.Echo "Exporting " & relativeBase
    doc.ExportAsFixedFormat freshPdf, 17
    doc.Close False
    fso.MoveFile freshPdf, finalPdf
    WScript.Echo "Completed " & relativeBase
Next
app.Quit
