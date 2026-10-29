' Read-only PDF export in a separate native Word instance; no study computation.
Option Explicit
Dim fso, taskRoot, app, doc, relativeBase, fileBase, targets
Set fso = CreateObject("Scripting.FileSystemObject")
taskRoot = fso.GetParentFolderName(fso.GetParentFolderName(WScript.ScriptFullName))
Set app = CreateObject("Word.Application")
If app.Documents.Count <> 0 Then Err.Raise 513, , "Expected an isolated Word instance"
app.Visible = False
app.DisplayAlerts = 0
targets = Array("manuscript_v1\submission_final_v2\Manuscript_BMC_Genomics_submission_v2", "manuscript_v1\submission_final_v2_zh\Manuscript_BMC_Genomics_submission_v2_zh", "manuscript_v1\submission_final_v2\Cover_letter_BMC_Genomics")
If WScript.Arguments.Count = 1 Then
    If WScript.Arguments(0) <> "--cover-only" Then Err.Raise 513, , "Unknown export selection"
    targets = Array("manuscript_v1\submission_final_v2\Cover_letter_BMC_Genomics")
End If
For Each relativeBase In targets
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
