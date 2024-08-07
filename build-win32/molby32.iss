[Setup]
AppName = Molby
AppVerName = Molby (v1.2beta1)
DefaultDirName = {pf}\Molby (32bit)
DefaultGroupName = Molby (32bit)
UninstallDisplayIcon = {app}\Molby.exe
OutputBaseFileName = SetupMolbyWin32

[Files]
Source: "build\release\Molby\Molby.exe"; DestDir: {app}
Source: "build\release\Molby\amber11\bin\*"; DestDir: {app}\amber11\bin
Source: "build\release\Molby\amber11\dat\antechamber\*"; DestDir: {app}\amber11\dat\antechamber
Source: "build\release\Molby\amber11\dat\leap\parm\*"; DestDir: {app}\amber11\dat\leap\parm
Source: "build\release\Molby\Scripts\*"; DestDir: {app}\Scripts
Source: "build\release\Molby\Scripts\mbsf\*"; DestDir: {app}\Scripts\mbsf
Source: "build\release\Molby\Scripts\mbsf\alicyclic\*"; DestDir: {app}\Scripts\mbsf\alicyclic
Source: "build\release\Molby\Scripts\mbsf\aromatic\*"; DestDir: {app}\Scripts\mbsf\aromatic
Source: "build\release\Molby\Scripts\mbsf\heterocyclic\*"; DestDir: {app}\Scripts\mbsf\heterocyclic
Source: "build\release\Molby\Scripts\mbsf\coordination\*"; DestDir: {app}\Scripts\mbsf\coordination
Source: "build\release\Molby\Scripts\mbsf\solvents\*"; DestDir: {app}\Scripts\mbsf\solvents
Source: "build\release\Molby\Scripts\mbsf\fragments\*"; DestDir: {app}\Scripts\mbsf\fragments
Source: "build\release\Molby\Scripts\lib\*"; DestDir: {app}\Scripts\lib
Source: "build\release\Molby\Scripts\basis_sets\*"; DestDir: {app}\Scripts\basis_sets
Source: "build\release\Molby\bitmaps\*"; DestDir: {app}\bitmaps
Source: "build\release\Molby\ortep3\*"; DestDir: {app}\ortep3

[Icons]
Name: "{group}\Molby"; Filename: "{app}\Molby.exe"
