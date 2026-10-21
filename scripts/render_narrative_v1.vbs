' Microsoft Word render only. Paths derive from this script, with no study-data inputs.
Option Explicit
Dim fso, taskRoot, taskOut, app, doc, stem, freshPdf, finalPdf
Set fso = CreateObject("Scripting.FileSystemObject")
taskRoot = fso.GetParentFolderName(fso.GetParentFolderName(WScript.ScriptFullName))
taskOut = fso.BuildPath(taskRoot, "manuscript_v1\author_review_v2")
Set app = CreateObject("Word.Application")
app.Visible = False
app.DisplayAlerts = 0
For Each stem In Array("Manuscript_BMC_Genomics_author_review_v2", "Supplementary_methods_author_review_v2")
    freshPdf = fso.BuildPath(taskOut, stem & ".vb_rendered.pdf")
    finalPdf = fso.BuildPath(taskOut, stem & ".pdf")
    If fso.FileExists(freshPdf) Then Err.Raise 513, , "Fresh render target already exists"
    WScript.Echo "Opening " & stem
    Set doc = app.Documents.Open(fso.BuildPath(taskOut, stem & ".docx"), False, True)
    WScript.Echo "Exporting " & stem
    doc.ExportAsFixedFormat freshPdf, 17
    doc.Close False
    Set doc = Nothing
    fso.CopyFile freshPdf, finalPdf, True
    fso.DeleteFile freshPdf, False
    WScript.Echo "Completed " & stem
Next
app.Quit
Set app = Nothing
