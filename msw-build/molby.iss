[Setup]
AppName = Molby
AppVerName = Molby (v0.5.0 build 20100121)
DefaultDirName = {pf}\Molby
DefaultGroupName = Molby
UninstallDisplayIcon = {app}\Molby.exe
OutputBaseFileName = SetupMolby

[Files]
Source: "build\Molby\Molby.exe"; DestDir: {app}
Source: "build\Molby\mingwm10.dll"; DestDir: {app}
Source: "build\Molby\Scripts\*"; DestDir: {app}\Scripts
Source: "build\Molby\Scripts\lib\*"; DestDir: {app}\Scripts\lib

[Icons]
Name: "{group}\Molby"; Filename: "{app}\Molby.exe"
