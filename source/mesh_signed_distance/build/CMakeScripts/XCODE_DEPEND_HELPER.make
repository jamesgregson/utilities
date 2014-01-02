# DO NOT EDIT
# This makefile makes sure all linkable targets are
# up-to-date with anything they link to
default:
	echo "Do not invoke directly"

# For each target create a dummy rule so the target does not have to exist
/usr/local/lib/libmpfr.dylib:
/usr/local/lib/libgmp.dylib:
/Users/jgregson/Code/ThirdParty/CGAL-4.1_build/lib/libCGAL.dylib:
/usr/local/lib/libboost_thread-mt.dylib:
/usr/local/lib/libboost_system-mt.dylib:


# Rules to remove targets that are older than anything to which they
# link.  This forces Xcode to relink the targets from scratch.  It
# does not seem to check these dependencies itself.
PostBuild.mesh_signed_distance.Debug:
/Users/jgregson/Code/GithubProjects/utilities/source/mesh_signed_distance/build/Debug/libmesh_signed_distance.dylib:\
	/usr/local/lib/libmpfr.dylib\
	/usr/local/lib/libgmp.dylib\
	/Users/jgregson/Code/ThirdParty/CGAL-4.1_build/lib/libCGAL.dylib\
	/usr/local/lib/libboost_thread-mt.dylib\
	/usr/local/lib/libboost_system-mt.dylib
	/bin/rm -f /Users/jgregson/Code/GithubProjects/utilities/source/mesh_signed_distance/build/Debug/libmesh_signed_distance.dylib


PostBuild.mesh_signed_distance.Release:
/Users/jgregson/Code/GithubProjects/utilities/libs/libmesh_signed_distance.dylib:\
	/usr/local/lib/libmpfr.dylib\
	/usr/local/lib/libgmp.dylib\
	/Users/jgregson/Code/ThirdParty/CGAL-4.1_build/lib/libCGAL.dylib\
	/usr/local/lib/libboost_thread-mt.dylib\
	/usr/local/lib/libboost_system-mt.dylib
	/bin/rm -f /Users/jgregson/Code/GithubProjects/utilities/libs/libmesh_signed_distance.dylib


PostBuild.mesh_signed_distance.MinSizeRel:
/Users/jgregson/Code/GithubProjects/utilities/source/mesh_signed_distance/build/MinSizeRel/libmesh_signed_distance.dylib:\
	/usr/local/lib/libmpfr.dylib\
	/usr/local/lib/libgmp.dylib\
	/Users/jgregson/Code/ThirdParty/CGAL-4.1_build/lib/libCGAL.dylib\
	/usr/local/lib/libboost_thread-mt.dylib\
	/usr/local/lib/libboost_system-mt.dylib
	/bin/rm -f /Users/jgregson/Code/GithubProjects/utilities/source/mesh_signed_distance/build/MinSizeRel/libmesh_signed_distance.dylib


PostBuild.mesh_signed_distance.RelWithDebInfo:
/Users/jgregson/Code/GithubProjects/utilities/source/mesh_signed_distance/build/RelWithDebInfo/libmesh_signed_distance.dylib:\
	/usr/local/lib/libmpfr.dylib\
	/usr/local/lib/libgmp.dylib\
	/Users/jgregson/Code/ThirdParty/CGAL-4.1_build/lib/libCGAL.dylib\
	/usr/local/lib/libboost_thread-mt.dylib\
	/usr/local/lib/libboost_system-mt.dylib
	/bin/rm -f /Users/jgregson/Code/GithubProjects/utilities/source/mesh_signed_distance/build/RelWithDebInfo/libmesh_signed_distance.dylib


