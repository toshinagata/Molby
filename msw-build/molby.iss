[Setup]
AppName = Molby
AppVerName = Molby (v0.6.4)
DefaultDirName = {pf}\Molby
DefaultGroupName = Molby
UninstallDisplayIcon = {app}\Molby.exe
OutputBaseFileName = SetupMolby

[Files]
Source: "build\Molby\Molby.exe"; DestDir: {app}
Source: "build\Molby\mingwm10.dll"; DestDir: {app}
Source: "build\Molby\amber11\bin\*"; DestDir: {app}\amber11\bin
Source: "build\Molby\amber11\dat\antechamber\*"; DestDir: {app}\amber11\dat\antechamber
Source: "build\Molby\amber11\dat\leap\parm\*"; DestDir: {app}\amber11\dat\leap\parm
Source: "build\Molby\Scripts\*"; DestDir: {app}\Scripts
Source: "build\Molby\Scripts\mbsf\*"; DestDir: {app}\Scripts\mbsf
Source: "build\Molby\Scripts\mbsf\alicyclic\*"; DestDir: {app}\Scripts\mbsf\alicyclic
Source: "build\Molby\Scripts\mbsf\aromatic\*"; DestDir: {app}\Scripts\mbsf\aromatic
Source: "build\Molby\Scripts\mbsf\heterocyclic\*"; DestDir: {app}\Scripts\mbsf\heterocyclic
Source: "build\Molby\Scripts\mbsf\solvents\*"; DestDir: {app}\Scripts\mbsf\solvents
Source: "build\Molby\Scripts\lib\*"; DestDir: {app}\Scripts\lib

[Icons]
Name: "{group}\Molby"; Filename: "{app}\Molby.exe"
