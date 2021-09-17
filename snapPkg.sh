VERSION=2.0.4
VERSION_=`echo -n $VERSION | tr '.' '_'`
export VERSION
PLATFORM=bionic
#svn copy https://svn.met.no/snap/trunk \
#     https://svn.met.no/snap/tags/$VERSION -m "Release $VERSION" || exit 1
#git checkout -b "version${VERSION_}"
#git commit -a -m "Version $VERSION"
#git push origin "version${VERSION_}"
#git checkout master
#git merge version${VERSION_}
cd src
make dist
cd ..
tar xvfz snap-${VERSION}.tgz
mv snap-${VERSION}.tgz bsnap_${VERSION}.orig.tar.gz
cd snap-${VERSION}
cp -r ../debian.${PLATFORM} .
mv debian.${PLATFORM} debian
dch -v ${VERSION}-1 "new release"
dch -r ""
debuild --preserve-envvar VERSION -us -uc -sa


