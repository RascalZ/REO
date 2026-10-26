' Read-only Word export of the Chinese reading translation; no scientific input is opened.
Option Explicit
Dim fso, taskRoot, taskOut, app, doc, stem, freshPdf, finalPdf
Set fso = CreateObject("Scripting.FileSystemObject")
taskRoot = fso.GetParentFolderName(fso.GetParentFolderName(WScript.ScriptFullName))
taskOut = fso.BuildPath(taskRoot, "manuscript_v1\author_review_v2_zh")
stem = "Manuscript_BMC_Genomics_author_review_v2_zh"
freshPdf = fso.BuildPath(taskOut, stem & ".vb_rendered.pdf")
finalPdf = fso.BuildPath(taskOut, stem & ".pdf")
If fso.FileExists(freshPdf) Or fso.FileExists(finalPdf) Then Err.Raise 513, , "Render target already exists; preserve prior reading copies"
Set app = CreateObject("Word.Application")
app.Visible = False
app.DisplayAlerts = 0
WScript.Echo "Opening Chinese reading Word"
Set doc = app.Documents.Open(fso.BuildPath(taskOut, stem & ".docx"), False, True)
doc.Repaginate
doc.Fields.Update
WScript.Echo "Exporting Chinese reading PDF"
doc.ExportAsFixedFormat freshPdf, 17
doc.Close False
app.Quit
fso.MoveFile freshPdf, finalPdf
WScript.Echo "Completed Chinese reading PDF"
