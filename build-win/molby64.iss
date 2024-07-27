;  Create an installer that can install both 32bit and 64bit versions
[Setup]
AppName = Molby
AppVerName = Molby (v1.2beta1)
DefaultDirName = {pf}\Molby
DefaultGroupName = Molby
UninstallDisplayIcon = {app}\Molby.exe
OutputBaseFileName = SetupMolbyWin
ArchitecturesInstallIn64BitMode = x64

[Files]
; Install 64-bit version if running in 64-bit mode, 32-bit version otherwise.
; x86_64 files
Source: "build\release\Molby\Molby.exe"; DestDir: {app}; Check: Is64BitInstallMode
Source: "build\release\Molby\amber11\bin\*"; DestDir: {app}\amber11\bin
Source: "build\release\Molby\ortep3\*"; DestDir: {app}\ortep3
; i386 files (the first one should be marked 'solidbreak')
Source: "..\build-win32\build\release\Molby\Molby.exe"; DestDir: {app}; Check: not Is64BitInstallMode; Flags: solidbreak
Source: "..\build-win32\build\release\Molby\amber11\bin\*"; DestDir: {app}\amber11\bin; Check: not Is64BitInstallMode
Source: "..\build-win32\build\release\Molby\ortep3\*"; DestDir: {app}\ortep3; Check: not Is64BitInstallMode
; Common files here (the first one should be marked 'solidbreak')
Source: "build\release\Molby\amber11\dat\antechamber\*"; DestDir: {app}\amber11\dat\antechamber; Flags: solidbreak
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
Source: "build\release\Molby\MolbyDoc\en\*"; DestDir: {app}\MolbyDoc\en
Source: "build\release\Molby\MolbyDoc\en\molby_rb\*"; DestDir: {app}\MolbyDoc\en\molby_rb
Source: "build\release\Molby\MolbyDoc\ja\*"; DestDir: {app}\MolbyDoc\ja
Source: "build\release\Molby\MolbyDoc\ja\molby_rb\*"; DestDir: {app}\MolbyDoc\ja\molby_rb
Source: "build\release\Molby\MolbyDoc\etc\*"; DestDir: {app}\MolbyDoc\etc
[Icons]
Name: "{group}\Molby"; Filename: "{app}\Molby.exe"
