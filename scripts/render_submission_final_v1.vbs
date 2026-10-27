' Native Word read-only export in an isolated instance; never close user files.
Option Explicit
Dim fso, taskRoot, app, doc, relativeBase, fileBase
Set fso = CreateObject("Scripting.FileSystemObject")
taskRoot = fso.GetParentFolderName(fso.GetParentFolderName(WScript.ScriptFullName))
Set app = CreateObject("Word.Application")
If app.Documents.Count <> 0 Then Err.Raise 513, , "Expected an isolated Word instance"
app.Visible = False
app.DisplayAlerts = 0
For Each relativeBase In Array("manuscript_v1\submission_final_v1\Manuscript_BMC_Genomics_submission", "manuscript_v1\submission_final_v1\Cover_letter_BMC_Genomics")
    fileBase = fso.BuildPath(taskRoot, relativeBase)
    If fso.FileExists(fileBase & ".pdf") Then Err.Raise 513, , "Preserve previous PDF"
    WScript.Echo "Opening " & relativeBase
    Set doc = app.Documents.Open(fileBase & ".docx", False, True)
    doc.Repaginate
    doc.Fields.Update
    doc.ExportAsFixedFormat fileBase & ".pdf", 17
    doc.Close False
    WScript.Echo "Exported " & relativeBase
Next
app.Quit
