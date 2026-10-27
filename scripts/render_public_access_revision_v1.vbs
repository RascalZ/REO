' Native Word read-only PDF export; never attaches to or closes the user's documents.
Option Explicit
Dim fso, taskRoot, app, doc, relativeBase, fileBase, freshPdf, finalPdf
Set fso = CreateObject("Scripting.FileSystemObject")
taskRoot = fso.GetParentFolderName(fso.GetParentFolderName(WScript.ScriptFullName))
Set app = CreateObject("Word.Application")
If app.Documents.Count <> 0 Then Err.Raise 513, , "Expected an isolated Word automation instance"
app.Visible = False
app.DisplayAlerts = 0
For Each relativeBase In Array("manuscript_v1\author_review_v5\Manuscript_BMC_Genomics_author_review_v5", "manuscript_v1\author_review_v5_zh\Manuscript_BMC_Genomics_author_review_v5_zh")
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
