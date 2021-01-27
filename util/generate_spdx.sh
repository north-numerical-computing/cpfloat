#! /bin/zsh

VERSION='0.1.0'

# Compute package hash using spdx-tools.
CURL=curl
INCDIR=./include/
SPDX_TOOLS_JAR='spdx-tools.jar'
SPDX_TOOLS_ULR='https://github.com/spdx/tools/releases/download/v2.2.4/spdx-tools-2.2.4-jar-with-dependencies.jar'
if [ ! -f $INCDIR$SPDX_TOOLS_JAR ]; then \
    $CURL -sL -o $INCDIR$SPDX_TOOLS_JAR $SPDX_TOOLS_ULR
fi

OUTPUT=`java -jar $INCDIR$SPDX_TOOLS_JAR \
      GenerateVerificationCode . ".*\.spdx|.*/include/.*|.*\.sh" \
      | awk -F ' ' 'NR==1{print $4}'`


# Add document and package information.
echo "##
## Document Creation Information
##

SPDXVersion: SPDX-2.2
DataLicense: CC0-1.0
SPDXID: SPDXRef-DOCUMENT
DocumentName: cpfloat-$VERSION
DocumentNamespace: https://raw.githubusercontent.com/mfasi/cpfloat/master/license.spdx
Creator: Person: Massimiliano Fasi (massimiliano.fasi@oru.se)
Creator: Person: Mantas Mikaitis (mantas.mikaitis@manchester.ac.uk)
Created: `date -u +%Y-%m-%dT%H:%M:%SZ`



##
## Package Information
##

PackageName: cpfloat
SPDXID: SPDXRef-1
PackageVersion: $VERSION
PackageDownloadLocation: git://github.com/mfasi/cpfloat
PackageVerificationCode: $OUTPUT (excludes: ./license.spdx)
PackageHomePage: https://github.com/mfasi/cpfloat
PackageLicenseConcluded: LGPL-2.1-or-later
PackageLicenseInfoFromFiles: LGPL-2.1-or-later
PackageLicenseDeclared: LGPL-2.1-or-later
PackageCopyrightText: <text>Copyright 2020 Massimiliano Fasi and Mantas Mikaitis</text>
PackageSummary:<text>Custom Precision Floating-point numbers.</text>



##
## File Information
##"

# Add file information.
counter=1
for file in `find .`; do
    if [[ ! -d $file && $file != (./.git*|./include/*|*.spdx) ]]; then
        echo ""
        echo "FileName: $file"
        echo "SPDXID: SPDXRef-1-$counter"
        counter=$((counter+1))
        case $file in
            *.sh|.git*|*/include/*|license.spdx)
            # Ignore:
            # * housekeeping scripts;
            # * git files;
            # * third-party files;
            # * license.spdx file.
            ;;
            *doc*|Doxyfile|cpfloat.m)
                echo "FileType: DOCUMENTATION"
                ;;
            Makefile|*.c|*.h|*.ts|*.m|*.cpp)
                echo "FileType: SOURCE"
                ;;
            *.md|*.txt)
                echo "FileType: TEXT"
                ;;
            *.spdx)
                echo "FileType: SPDX"
                ;;
            *)
                echo "FileType: OTHER"
        esac
        echo "FileChecksum: SHA1: `shasum -a 1 $file | \
                                      awk -F ' ' '{print $1}'`"
        echo "FileChecksum: MD5: `md5 $file | awk -F ' ' '{print $4}'`"
        echo "LicenseConcluded: LGPL-2.1-or-later"
        LICENSE=`grep "SPDX-License-Identifier" $file | \
                      awk -F ' ' '{printf $3}'`
        if [[ $LICENSE = "" ]]; then
            LICENSE=NONE
        fi
        echo "LicenseInfoInFile: $LICENSE"
    fi
done
