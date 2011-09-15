[Setup]
AppName = Molby
AppVerName = Molby (v0.5.6.2)
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
Source: "build\Molby\Scripts\lib\*"; DestDir: {app}\Scripts\lib

[Icons]
Name: "{group}\Molby"; Filename: "{app}\Molby.exe"
