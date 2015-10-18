#!/bin/bash

mkdir ./xmolview.app
mkdir ./xmolview.app/Contents
mkdir ./xmolview.app/Contents/Resources
mkdir ./xmolview.app/Contents/MacOS
echo APPLnone > ./xmolview.app/Contents/PkgInfo
cp bin/xmolview ./xmolview.app/Contents/MacOS
chmod 755 ./xmolview.app/Contents/MacOS/xmolview
cat << EOF > ./xmolview.app/Contents/info.plist
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist SYSTEM "file://localhost/System/Library/DTDs/PropertyList.dtd">
<plist version="0.9">
<dict>
        <key>CFBundleName</key>
        <string>xmolview</string>
        <key>CFBundlePackageType</key>
        <string>APPL</string>
        <key>CFBundleVersion</key>
        <string>59</string>
        <key>CFBundleShortVersionString</key>
        <string>1.1</string>
        <key>CFBundleIconFile</key>
        <string>icon.icns</string>
        <key>CFBundleSignature</key>
        <string>none</string>
</dict>
</plist>
EOF

