/*
 * InterfaceI23dSFM.cpp
 *
 * Copyright (c) 2014-2015 I23D
 *
 *
 *      Pierre MOULON <p.moulon@foxel.ch>
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * Additional Terms:
 *
 *      You are required to preserve legal notices and author attributions in
 *      that material or in the Appropriate Legal Notices displayed by works
 *      containing it.
 */

#include "../../libs/MVS/Common.h"
#include "../../libs/MVS/Scene.h"
#include <boost/program_options.hpp>
#include <fstream>
#include <vector>
#include <set>
#include "MVSAllSteps.h"

#ifdef _USE_I23DSFM
#undef D2R
#undef R2D
#include <i23dSFM/sfm/sfm_data.hpp>
#include <i23dSFM/sfm/sfm_data_io.hpp>
#include <i23dSFM/image/image.hpp>
#endif

std::string sfm_data_path;
std::string output_dir;
std::string working_dir;
std::string I23D_BIN = "/edata/lyb/i23dMVS2_build/bin/";
// D E F I N E S ///////////////////////////////////////////////////

#define APPNAME _T("MVSAllSteps")
#define MVS_EXT _T(".mvs")
#define MVG_EXT _T(".baf")
#define MVG2_EXT _T(".json")

using namespace MVS;
// S T R U C T S ///////////////////////////////////////////////////

namespace i23D {
namespace MVS_IO {

typedef REAL RealT;
typedef Eigen::Matrix<RealT,3,3,Eigen::RowMajor> Mat33;
typedef Eigen::Matrix<RealT,3,1> Vec3;

// Structure to model the pinhole camera projection model
struct Camera
{
	Mat33 K; // camera's normalized intrinsics matrix
};
typedef std::vector<Camera> vec_Camera;

// structure describing a pose along the trajectory of a platform
struct Pose {
	Mat33 R;  // pose's rotation matrix
	Vec3 C;   // pose's translation vector
};
typedef std::vector<Pose> vec_Pose;

// structure describing an image
struct Image {
	uint32_t id_camera; // ID of the associated camera on the associated platform
	uint32_t id_pose;   // ID of the pose of the associated platform
	std::string name;   // image file name
};
typedef std::vector<Image> vec_Image;

// structure describing a 3D point
struct Vertex {

	typedef std::vector<uint32_t> vec_View;

	Vec3 X; // 3D point position
	vec_View views; // view visibility for this 3D feature
};
typedef std::vector<Vertex> vec_Vertex;

struct SfM_Scene
{
	vec_Pose poses;       // array of poses
	vec_Camera cameras;   // array of cameras
	vec_Image images;     // array of images
	vec_Vertex vertices;  // array of reconstructed 3D points
};


bool ImportScene(const std::string& sList_filename, const std::string& sBaf_filename, SfM_Scene& sceneBAF)
{
	LOG_OUT() << "Reading:\n"
		<< sList_filename << "\n"
		<< sBaf_filename << std::endl;

	// Read view list file (view filename, id_intrinsic, id_pose)
	// Must be read first, since it allow to establish the link between the ViewId and the camera/poses ids.
	std::map< std::pair<uint32_t, uint32_t>, uint32_t > map_cam_pose_toViewId;
	{
		std::ifstream file(sList_filename.c_str());
		if (!file.good()) {
			VERBOSE("error: unable to open file '%s'", sList_filename.c_str());
			return false;
		}
		Image image;
		uint32_t count = 0;
		while (file >> image.name >> image.id_camera >> image.id_pose) {
			sceneBAF.images.push_back(image);
			map_cam_pose_toViewId[std::make_pair(image.id_camera, image.id_pose)] = count++;
			LOG_OUT() << image.name << ' ' << image.id_camera << ' ' << image.id_pose << std::endl;
		}
	}

	// Read BAF file
	{
		std::ifstream file(sBaf_filename.c_str());
		if (!file.good()) {
			VERBOSE("error: unable to open file '%s'", sBaf_filename.c_str());
			return false;
		}

		uint32_t num_intrinsics, num_poses, num_points;

		// Read header
		file >> num_intrinsics;
		file >> num_poses;
		file >> num_points;

		LOG_OUT() << "Reading BAF file with:\n"
			<< " num_intrinsics: " << num_intrinsics << "\n"
			<< " num_poses: " << num_poses << "\n"
			<< " num_points: " << num_points << "\n";

		// Read the intrinsics (only support reading Pinhole Radial 3).
		{
			for (uint32_t i = 0; i < num_intrinsics; ++i) {
				double focal, ppx, ppy, k1, k2, k3;
				file >> focal >> ppx >> ppy >> k1 >> k2 >> k3;
				Camera cam;
				cam.K <<
					focal, 0, ppx,
					0, focal, ppy,
					0, 0, 1;
				LOG_OUT() << "\n" << cam.K << std::endl;
				sceneBAF.cameras.push_back(cam);
			}
		}

		// Read poses
		{
			for (uint32_t i = 0; i < num_poses; ++i) {
				Pose pose;
				for (int r = 0; r < 3; ++r) {
					for (int c = 0; c < 3; ++c) {
						file >> pose.R(r,c);
					}
				}
				file >> pose.C[0] >> pose.C[1] >> pose.C[2];
				#ifndef _RELEASE
				LOG_OUT() << "\n" << pose.R << "\n\n" << pose.C.transpose() << std::endl;
				#endif
				sceneBAF.poses.push_back(pose);
			}
		}

		// Read structure and visibility
		{
			#ifdef _RELEASE
			Util::Progress progress(_T("Processed points"), num_points);
			#endif
			for (uint32_t i = 0; i < num_points; ++i) {
				Vertex vertex;
				file >> vertex.X[0] >> vertex.X[1] >> vertex.X[2];
				uint32_t num_observations_for_point = 0;
				file >> num_observations_for_point;
				for (uint32_t j = 0; j < num_observations_for_point; ++j) {
					uint32_t id_intrinsics, id_pose;
					double x, y;
					file >> id_intrinsics >> id_pose >> x >> y;
					#ifndef _RELEASE
					LOG_OUT() << "observation:"
						<< " " <<  id_intrinsics
						<< " " <<  id_pose
						<< " " << x << " " << y << std::endl;
					#endif
					const auto itIntrPose(map_cam_pose_toViewId.find(std::make_pair(id_intrinsics, id_pose)));
					if (itIntrPose == map_cam_pose_toViewId.end()) {
						LOG_OUT() << "error: intrinsics-pose pair not existing" << std::endl;
						continue;
					}
					const uint32_t id_view(itIntrPose->second);
					vertex.views.push_back(id_view);
				}
				sceneBAF.vertices.push_back(vertex);
				#ifdef _RELEASE
				progress.display(i);
				#endif
			}
			#ifdef _RELEASE
			progress.close();
			#endif
		}
	}
	return true;
}

bool ExportScene(const std::string& sList_filename, const std::string& sBaf_filename, const SfM_Scene& sceneBAF)
{
	LOG_OUT() << "Writing:\n"
		<< sList_filename << "\n"
		<< sBaf_filename << std::endl;

	// Write view list file (view filename, id_intrinsic, id_pose)
	{
		std::ofstream file(sList_filename.c_str());
		if (!file.good()) {
			VERBOSE("error: unable to open file '%s'", sList_filename.c_str());
			return false;
		}
		for (uint32_t i=0; i<sceneBAF.images.size(); ++i) {
			const Image& image = sceneBAF.images[i];
			file << image.name << ' ' << image.id_camera << ' ' << image.id_pose << std::endl;
			LOG_OUT() << image.name << ' ' << image.id_camera << ' ' << image.id_pose << std::endl;
		}
	}

	// Write BAF file
	{
		std::ofstream file(sBaf_filename.c_str());
		if (!file.good()) {
			VERBOSE("error: unable to open file '%s'", sBaf_filename.c_str());
			return false;
		}

		const uint32_t num_intrinsics = (uint32_t)sceneBAF.cameras.size();
		const uint32_t num_poses = (uint32_t)sceneBAF.poses.size();
		const uint32_t num_points = (uint32_t)sceneBAF.vertices.size();

		LOG_OUT() << "Writing BAF file with:\n"
			<< " num_intrinsics: " << num_intrinsics << "\n"
			<< " num_poses: " << num_poses << "\n"
			<< " num_points: " << num_points << "\n";

		// Write header
		file << num_intrinsics << std::endl;
		file << num_poses << std::endl;
		file << num_points << std::endl;

		// Write the intrinsics (only support writing Pinhole Radial 3).
		{
			for (uint32_t i = 0; i < num_intrinsics; ++i) {
				const Camera& cam = sceneBAF.cameras[i];
				file << cam.K(0,0) << ' ' << cam.K(0,2) << ' ' << cam.K(1,2) << ' ' << 0 << ' ' << 0 << ' ' << 0 << std::endl;
				LOG_OUT() << "\n" << cam.K << std::endl;
			}
		}

		// Write poses
		{
			for (uint32_t i = 0; i < num_poses; ++i) {
				const Pose& pose = sceneBAF.poses[i];
				for (int r = 0; r < 3; ++r) {
					for (int c = 0; c < 3; ++c) {
						file << pose.R(r,c) << ' ';
					}
				}
				file << pose.C[0] << ' ' << pose.C[1] << ' ' << pose.C[2] << std::endl;
				#ifndef _RELEASE
				LOG_OUT() << "\n" << pose.R << "\n\n" << pose.C.transpose() << std::endl;
				#endif
			}
		}

		// Write structure and visibility
		{
			#ifdef _RELEASE
			Util::Progress progress(_T("Processed points"), num_points);
			#endif
			for (uint32_t i = 0; i < num_points; ++i) {
				const Vertex& vertex = sceneBAF.vertices[i];
				file << vertex.X[0] << ' ' << vertex.X[1] << ' ' << vertex.X[2] << std::endl;
				const uint32_t num_observations_for_point = (uint32_t)vertex.views.size();
				file << num_observations_for_point << std::endl;
				for (uint32_t j = 0; j < num_observations_for_point; ++j) {
					const uint32_t id_view = vertex.views[j];
					const Image& image = sceneBAF.images[id_view];
					file << image.id_camera << ' ' << image.id_pose << ' ' << 0 << ' ' << 0 << std::endl;
					#ifndef _RELEASE
					LOG_OUT() << "observation:"
						<< " " <<  image.id_camera
						<< " " <<  image.id_pose << std::endl;
					#endif
				}
				#ifdef _RELEASE
				progress.display(i);
				#endif
			}
			#ifdef _RELEASE
			progress.close();
			#endif
		}
	}
	return true;
}
} // MVS_IO
} // i23D


namespace OPT {
#ifdef _USE_I23DSFM
bool bI23dSFMjson; // new import format
#endif
bool bI23D2I23dSFM; // conversion direction
bool bNormalizeIntrinsics;
String strListFileName;
String strInputFileName;
String strOutputFileName;
String strOutputImageFolder;
unsigned nArchiveType;
int nProcessPriority;
unsigned nMaxThreads;
String strConfigFileName;
boost::program_options::variables_map vm;

String strMeshFileName;
String strDenseConfigFileName;
unsigned nResolutionLevel;
unsigned nMinResolution;
unsigned nEstimateColors;
unsigned nEstimateNormals;

bool bMeshExport;
float fDistInsert;
bool bUseFreeSpaceSupport;
float fDecimateMesh;
float fRemoveSpurious;
bool bRemoveSpikes;
unsigned nCloseHoles;
unsigned nSmoothMesh;

unsigned nMaxViews;
unsigned nEnsureEdgeSize;
unsigned nScales;
float fScaleStep;
unsigned nReduceMemory;
unsigned nAlternatePair;
float fRegularityWeight;
float fRatioRigidityElasticity;
unsigned nMaxFaceArea;
float fPlanarVertexRatio;
float fGradientStep;
#ifdef _USE_CUDA
bool bUseCUDA;
#endif

float fOutlierThreshold;
float fRatioDataSmoothness;
bool bGlobalSeamLeveling;
bool bLocalSeamLeveling;
unsigned nTextureSizeMultiple;
unsigned nRectPackingHeuristic;
uint32_t nColEmpty;

String strViwoConfigFileName;
std::vector<String> inputFileNames;
std::vector<String> transFileNames;
} // namespace OPT

// initialize and parse the command line parameters
bool Initialize_InterfaceI23dSFM()
{
/*
	size_t argc=7;
	const char * argv[8];
	//for (int i=0;i<7;i++) argv[i]=new char[100];
	argv[0]="InterfaceI23dSFM";
	argv[1]="-i";
	argv[2]=sfm_data_path.c_str();
	argv[3]="-o";
	argv[4]=(output_dir+"/scene.mvs").c_str();
	argv[5]="-w";
	argv[6]=working_dir.c_str();
	argv[argc]=NULL;
*/
	OPT::strInputFileName=sfm_data_path;
	OPT::strOutputFileName=output_dir+"/scene.mvs";
	WORKING_FOLDER=working_dir;
	OPT::strConfigFileName="InterfaceI23dSFM.cfg";
	OPT::nArchiveType=2;
	OPT::nProcessPriority=-1;
	OPT::nMaxThreads=0;
	g_nVerbosityLevel=2;
	OPT::strOutputImageFolder="undistorted_images";
	OPT::bNormalizeIntrinsics=true;
	// initialize log and console
	OPEN_LOG();
	OPEN_LOGCONSOLE();
	//for (int i=0;i<7;i++) std::cout<<argv[i]<<std::endl;
	// group of options allowed only on command line
	/*
	boost::program_options::options_description generic("Generic options");
	generic.add_options()
		("help,h", "produce this help message")
		("working-folder,w", boost::program_options::value<std::string>(&WORKING_FOLDER), "working directory (default current directory)")
		("config-file,c", boost::program_options::value<std::string>(&OPT::strConfigFileName)->default_value(APPNAME _T(".cfg")), "file name containing program options")
		("archive-type", boost::program_options::value<unsigned>(&OPT::nArchiveType)->default_value(2), "project archive type: 0-text, 1-binary, 2-compressed binary")
		("process-priority", boost::program_options::value<int>(&OPT::nProcessPriority)->default_value(-1), "process priority (below normal by default)")
		("max-threads", boost::program_options::value<unsigned>(&OPT::nMaxThreads)->default_value(0), "maximum number of threads (0 for using all available cores)")
		#if TD_VERBOSE != TD_VERBOSE_OFF
		("verbosity,v", boost::program_options::value<int>(&g_nVerbosityLevel)->default_value(
			#if TD_VERBOSE == TD_VERBOSE_DEBUG
			3
			#else
			2
			#endif
			), "verbosity level")
		#endif
		;
	// group of options allowed both on command line and in config file
	boost::program_options::options_description config("Main options");
	config.add_options()
		("images-list-file,l", boost::program_options::value<std::string>(&OPT::strListFileName), "input filename containing image list")
		("input-file,i", boost::program_options::value<std::string>(&OPT::strInputFileName), "input filename containing camera poses and image list")
		("output-file,o", boost::program_options::value<std::string>(&OPT::strOutputFileName), "output filename for storing the mesh")
		("output-image-folder", boost::program_options::value<std::string>(&OPT::strOutputImageFolder)->default_value("undistorted_images"), "output folder to store undistorted images")
		("normalize,f", boost::program_options::value<bool>(&OPT::bNormalizeIntrinsics)->default_value(true), "normalize intrinsics while exporting to I23D format")
		;

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(config);

	boost::program_options::options_description config_file_options;
	config_file_options.add(config);

	boost::program_options::positional_options_description p;
	p.add("input-file", -1);

	try {
		// parse command line options
		boost::program_options::store(boost::program_options::command_line_parser((int)argc, argv).options(cmdline_options).positional(p).run(), OPT::vm);
		boost::program_options::notify(OPT::vm);
		INIT_WORKING_FOLDER;
		// parse configuration file
		std::ifstream ifs(MAKE_PATH_SAFE(OPT::strConfigFileName));
		//fault
		if (ifs) {
			boost::program_options::store(parse_config_file(ifs, config_file_options), OPT::vm);
			boost::program_options::notify(OPT::vm);
		}
	}
	catch (const std::exception& e) {
		LOG(e.what());
		return false;
	}
	*/
	//std::cout<<"88888        "<<OPT::strInputFileName<<std::endl;
    //std::cout<<OPT::strOutputFileName<<std::endl;
	// initialize the log file
	INIT_WORKING_FOLDER;
	OPEN_LOGFILE(MAKE_PATH(APPNAME _T("-")+Util::getUniqueName(0)+_T(".log")).c_str());

	// print application details: version and command line
	Util::LogBuild();
	//std::cout<<Util::CommandLineToString(argc, argv).c_str()<<std::endl;
	//for (int i=0;i<7;i++) std::cout<<argv[i]<<std::endl;
	//LOG(_T("Command line:%s"), Util::CommandLineToString(argc, argv).c_str());

	// validate input
	Util::ensureValidPath(OPT::strListFileName);
	Util::ensureUnifySlash(OPT::strListFileName);
	Util::ensureValidPath(OPT::strInputFileName);
	Util::ensureUnifySlash(OPT::strInputFileName);
	Util::ensureUnifySlash(OPT::strOutputImageFolder);
	Util::ensureDirectorySlash(OPT::strOutputImageFolder);
	const String strInputFileNameExt(Util::getFileExt(OPT::strInputFileName).ToLower());
	OPT::bI23D2I23dSFM = (strInputFileNameExt == MVS_EXT);
	#ifdef _USE_I23DSFM
	OPT::bI23dSFMjson = (strInputFileNameExt == MVG2_EXT);
	const bool bInvalidCommand(OPT::strInputFileName.IsEmpty() || (OPT::strListFileName.IsEmpty() && !OPT::bI23dSFMjson && !OPT::bI23D2I23dSFM));
	#else
	const bool bInvalidCommand(OPT::strInputFileName.IsEmpty() || (OPT::strListFileName.IsEmpty() && !OPT::bI23D2I23dSFM));
	#endif
	//if (OPT::vm.count("help") || bInvalidCommand) {
		//boost::program_options::options_description visible("Available options");
		//visible.add(generic).add(config);
		//GET_LOG() << visible;
	//}
	if (bInvalidCommand)
		return false;

	// initialize optional options
	Util::ensureValidPath(OPT::strOutputFileName);
	Util::ensureUnifySlash(OPT::strOutputFileName);
	if (OPT::bI23D2I23dSFM) {
		if (OPT::strOutputFileName.IsEmpty())
			OPT::strOutputFileName = Util::getFullFileName(OPT::strInputFileName);
	} else {
		if (OPT::strOutputFileName.IsEmpty())
			OPT::strOutputFileName = Util::getFullFileName(OPT::strInputFileName) + MVS_EXT;
	}

	// initialize global options
	Process::setCurrentProcessPriority((Process::Priority)OPT::nProcessPriority);
	#ifdef _USE_OPENMP
	if (OPT::nMaxThreads != 0)
		omp_set_num_threads(OPT::nMaxThreads);
	#endif

	#ifdef _USE_BREAKPAD
	// start memory dumper
	MiniDumper::Create(APPNAME, WORKING_FOLDER);
	#endif
	return true;
}
// initialize and parse the command line parameters
bool Initialize_DensifyPointCloud()
{
	/*
	size_t argc=7;
	const char * argv[8];
	argv[0]="MVSAllSteps";
	argv[1]="-i";
	argv[2]=(output_dir+"/scene.mvs").c_str();
	argv[3]="-w";
	argv[4]=working_dir.c_str();
	argv[5]="--estimate-normals";
	argv[6]="1";
	//WORKING_FOLDER=working_dir;
	//OPT::strInputFileName=output_dir+"/scene.mvs";
	argv[argc]=NULL;
	*/
	// initialize log and console
	OPEN_LOG();
	OPEN_LOGCONSOLE();
	
	WORKING_FOLDER=working_dir;
	OPT::strInputFileName=output_dir+"/scene.mvs";
	OPT::strConfigFileName="DensifyPointCloud.cfg";
	OPT::nArchiveType=2;
	OPT::nProcessPriority=2;
	OPT::nMaxThreads=0;
	g_nVerbosityLevel=1;
	OPT::nResolutionLevel=1;
	OPT::nMinResolution=640;
	OPT::nEstimateColors=1;
	OPT::nEstimateNormals=2;
	/*
	// group of options allowed only on command line
	boost::program_options::options_description generic("Generic options");
	generic.add_options()
		("help,h", "produce this help message")
		("working-folder,w", boost::program_options::value<std::string>(&WORKING_FOLDER), "working directory (default current directory)")
		("config-file,c", boost::program_options::value<std::string>(&OPT::strConfigFileName)->default_value(APPNAME _T(".cfg")), "file name containing program options")
		("archive-type", boost::program_options::value<unsigned>(&OPT::nArchiveType)->default_value(2), "project archive type: 0-text, 1-binary, 2-compressed binary")
		("process-priority", boost::program_options::value<int>(&OPT::nProcessPriority)->default_value(2), "process priority (below normal by default)")
		("max-threads", boost::program_options::value<unsigned>(&OPT::nMaxThreads)->default_value(0), "maximum number of threads (0 for using all available cores)")
		//#if TD_VERBOSE != TD_VERBOSE_OFF
		("verbosity,v", boost::program_options::value<int>(&g_nVerbosityLevel)->default_value(
			#if TD_VERBOSE == TD_VERBOSE_DEBUG
			3
			#else
			1
			#endif
			), "verbosity level")
		//#endif
		;

	// group of options allowed both on command line and in config file
	boost::program_options::options_description config("Refine options");
	config.add_options()
		("input-file,i", boost::program_options::value<std::string>(&OPT::strInputFileName), "input filename containing camera poses and image list")
		("output-file,o", boost::program_options::value<std::string>(&OPT::strOutputFileName), "output filename for storing the dense point-cloud")
		("resolution-level", boost::program_options::value<unsigned>(&OPT::nResolutionLevel)->default_value(1), "how many times to scale down the images before point cloud computation")
		("min-resolution", boost::program_options::value<unsigned>(&OPT::nMinResolution)->default_value(640), "do not scale images lower than this resolution")
		("estimate-colors", boost::program_options::value<unsigned>(&OPT::nEstimateColors)->default_value(1), "estimate the colors for the dense point-cloud")
		("estimate-normals", boost::program_options::value<unsigned>(&OPT::nEstimateNormals)->default_value(2), "estimate the normals for the dense point-cloud")
		;

	// hidden options, allowed both on command line and
	// in config file, but will not be shown to the user
	boost::program_options::options_description hidden("Hidden options");
	hidden.add_options()
		("dense-config-file", boost::program_options::value<std::string>(&OPT::strDenseConfigFileName), "optional configuration file for the densifier (overwritten by the command line options)")
		;

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(config).add(hidden);

	boost::program_options::options_description config_file_options;
	config_file_options.add(config).add(hidden);

	boost::program_options::positional_options_description p;
	p.add("input-file", -1);

	try {
		// parse command line options
		boost::program_options::store(boost::program_options::command_line_parser((int)argc, argv).options(cmdline_options).positional(p).run(), OPT::vm);
		boost::program_options::notify(OPT::vm);
		INIT_WORKING_FOLDER;
		// parse configuration file
		std::ifstream ifs(MAKE_PATH_SAFE(OPT::strConfigFileName));
		if (ifs) {
			boost::program_options::store(parse_config_file(ifs, config_file_options), OPT::vm);
			boost::program_options::notify(OPT::vm);
		}
	}
	catch (const std::exception& e) {
		LOG(e.what());
		return false;
	}
	*/
	// initialize the log file
	INIT_WORKING_FOLDER;

	OPEN_LOGFILE(MAKE_PATH(APPNAME _T("-")+Util::getUniqueName(0)+_T(".log")).c_str());

	// print application details: version and command line
	Util::LogBuild();
	//LOG(_T("Command line:%s"), Util::CommandLineToString(argc, argv).c_str());

	// validate input
	Util::ensureValidPath(OPT::strInputFileName);
	Util::ensureUnifySlash(OPT::strInputFileName);
	if (OPT::vm.count("help") || OPT::strInputFileName.IsEmpty()) {
		boost::program_options::options_description visible("Available options");
		//visible.add(generic).add(config);
		//GET_LOG() << visible;
	}
	if (OPT::strInputFileName.IsEmpty())
		return false;

	// initialize optional options
	Util::ensureValidPath(OPT::strOutputFileName);
	Util::ensureUnifySlash(OPT::strOutputFileName);
	if (OPT::strOutputFileName.IsEmpty())
		OPT::strOutputFileName = Util::getFullFileName(OPT::strInputFileName) + _T(".mvs");

	// init dense options
	OPTDENSE::init();
	const bool bValidConfig(OPTDENSE::oConfig.Load(OPT::strDenseConfigFileName));
	OPTDENSE::update();
	OPTDENSE::nResolutionLevel = OPT::nResolutionLevel;
	OPTDENSE::nMinResolution = OPT::nMinResolution;
	OPTDENSE::nEstimateColors = OPT::nEstimateColors;
	OPTDENSE::nEstimateNormals = OPT::nEstimateNormals;
	if (!bValidConfig)
		OPTDENSE::oConfig.Save(OPT::strDenseConfigFileName);

	// initialize global options
	Process::setCurrentProcessPriority((Process::Priority)OPT::nProcessPriority);
	#ifdef _USE_OPENMP
	if (OPT::nMaxThreads != 0)
		omp_set_num_threads(OPT::nMaxThreads);
	#endif

	#ifdef _USE_BREAKPAD
	// start memory dumper
	MiniDumper::Create(APPNAME, WORKING_FOLDER);
	#endif
	return true;
}
// initialize and parse the command line parameters
bool Initialize_ReconstructMesh()
{
	size_t argc=5;
	const char * argv[5];
	//for (int i=0;i<7;i++) argv[i]=new char[100];
	argv[0]=(I23D_BIN+"MVSAllSteps").c_str();
	argv[1]="-i";
	argv[2]=(output_dir+"/scene_dense.mvs").c_str();
	argv[3]="-w";
	argv[4]=working_dir.c_str();

	// initialize log and console
	OPEN_LOG();
	OPEN_LOGCONSOLE();

	// group of options allowed only on command line
	boost::program_options::options_description generic("Generic options");
	generic.add_options()
		("help,h", "produce this help message")
		("working-folder,w", boost::program_options::value<std::string>(&WORKING_FOLDER), "working directory (default current directory)")
		("config-file,c", boost::program_options::value<std::string>(&OPT::strConfigFileName)->default_value(APPNAME _T(".cfg")), "file name containing program options")
		("archive-type", boost::program_options::value<unsigned>(&OPT::nArchiveType)->default_value(2), "project archive type: 0-text, 1-binary, 2-compressed binary")
		("process-priority", boost::program_options::value<int>(&OPT::nProcessPriority)->default_value(-1), "process priority (below normal by default)")
		("max-threads", boost::program_options::value<unsigned>(&OPT::nMaxThreads)->default_value(0), "maximum number of threads (0 for using all available cores)")
		#if TD_VERBOSE != TD_VERBOSE_OFF
		("verbosity,v", boost::program_options::value<int>(&g_nVerbosityLevel)->default_value(
			#if TD_VERBOSE == TD_VERBOSE_DEBUG
			3
			#else
			2
			#endif
			), "verbosity level")
		#endif
		;

	// group of options allowed both on command line and in config file
	boost::program_options::options_description config_main("Reconstruct options");
	config_main.add_options()
		("input-file,i", boost::program_options::value<std::string>(&OPT::strInputFileName), "input filename containing camera poses and image list")
		("output-file,o", boost::program_options::value<std::string>(&OPT::strOutputFileName), "output filename for storing the mesh")
		("min-point-distance,d", boost::program_options::value<float>(&OPT::fDistInsert)->default_value(2.f), "minimum distance in pixels between the projection of two 3D points to consider them different while triangulating")
		("free-space-support,f", boost::program_options::value<bool>(&OPT::bUseFreeSpaceSupport)->default_value(true), "exploits the free-space support in order to reconstruct weakly-represented surfaces")
		;
	boost::program_options::options_description config_clean("Clean options");
	config_clean.add_options()
		("decimate", boost::program_options::value<float>(&OPT::fDecimateMesh)->default_value(1.f), "decimation factor in range (0..1] to be applied to the reconstructed surface (1 - disabled)")
		("remove-spurious", boost::program_options::value<float>(&OPT::fRemoveSpurious)->default_value(20.f), "spurious factor for removing faces with too long edges or isolated components (0 - disabled)")
		("remove-spikes", boost::program_options::value<bool>(&OPT::bRemoveSpikes)->default_value(true), "flag controlling the removal of spike faces")
		("close-holes", boost::program_options::value<unsigned>(&OPT::nCloseHoles)->default_value(30), "try to close small holes in the reconstructed surface (0 - disabled)")
		("smooth", boost::program_options::value<unsigned>(&OPT::nSmoothMesh)->default_value(2), "number of iterations to smooth the reconstructed surface (0 - disabled)")
		;

	// hidden options, allowed both on command line and
	// in config file, but will not be shown to the user
	boost::program_options::options_description hidden("Hidden options");
	hidden.add_options()
		("mesh-file", boost::program_options::value<std::string>(&OPT::strMeshFileName), "mesh file name to clean (skips the reconstruction step)")
		("mesh-export", boost::program_options::value<bool>(&OPT::bMeshExport)->default_value(false), "just export the mesh contained in loaded project")
		;

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(config_main).add(config_clean).add(hidden);

	boost::program_options::options_description config_file_options;
	config_file_options.add(config_main).add(config_clean).add(hidden);

	boost::program_options::positional_options_description p;
	p.add("input-file", -1);

	try {
		// parse command line options
		boost::program_options::store(boost::program_options::command_line_parser((int)argc, argv).options(cmdline_options).positional(p).run(), OPT::vm);
		boost::program_options::notify(OPT::vm);
		INIT_WORKING_FOLDER;
		// parse configuration file
		std::ifstream ifs(MAKE_PATH_SAFE(OPT::strConfigFileName));
		if (ifs) {
			boost::program_options::store(parse_config_file(ifs, config_file_options), OPT::vm);
			boost::program_options::notify(OPT::vm);
		}
	}
	catch (const std::exception& e) {
		LOG(e.what());
		return false;
	}

	// initialize the log file
	OPEN_LOGFILE(MAKE_PATH(APPNAME _T("-")+Util::getUniqueName(0)+_T(".log")).c_str());

	// print application details: version and command line
	Util::LogBuild();
	LOG(_T("Command line:%s"), Util::CommandLineToString(argc, argv).c_str());

	// validate input
	Util::ensureValidPath(OPT::strInputFileName);
	Util::ensureUnifySlash(OPT::strInputFileName);
	if (OPT::vm.count("help") || OPT::strInputFileName.IsEmpty()) {
		boost::program_options::options_description visible("Available options");
		visible.add(generic).add(config_main).add(config_clean);
		GET_LOG() << visible;
	}
	if (OPT::strInputFileName.IsEmpty())
		return false;

	// initialize optional options
	Util::ensureValidPath(OPT::strOutputFileName);
	Util::ensureUnifySlash(OPT::strOutputFileName);
	if (OPT::strOutputFileName.IsEmpty())
		OPT::strOutputFileName = Util::getFullFileName(OPT::strInputFileName) + _T(".mvs");

	// initialize global options
	Process::setCurrentProcessPriority((Process::Priority)OPT::nProcessPriority);
	#ifdef _USE_OPENMP
	if (OPT::nMaxThreads != 0)
		omp_set_num_threads(OPT::nMaxThreads);
	#endif

	#ifdef _USE_BREAKPAD
	// start memory dumper
	MiniDumper::Create(APPNAME, WORKING_FOLDER);
	#endif
	return true;
}
// initialize and parse the command line parameters
bool Initialize_RefineMesh()
{
	// initialize log and console
	OPEN_LOG();
	OPEN_LOGCONSOLE();
	size_t argc=9;
	const char * argv[9];
	//for (int i=0;i<9;i++) argv[i]=new char[100];
	argv[0]=(I23D_BIN+"MVSAllSteps").c_str();
	argv[1]="-i";
	argv[2]=(output_dir+"/scene_dense_mesh.mvs").c_str();
	argv[3]="-w";
	argv[4]=working_dir.c_str();
	argv[5]="--scales";
	argv[6]="2";
	argv[7]="--resolution-level";
	argv[8]="2";
	// group of options allowed only on command line
	boost::program_options::options_description generic("Generic options");
	generic.add_options()
		("help,h", "produce this help message")
		("working-folder,w", boost::program_options::value<std::string>(&WORKING_FOLDER), "working directory (default current directory)")
		("config-file,c", boost::program_options::value<std::string>(&OPT::strConfigFileName)->default_value(APPNAME _T(".cfg")), "file name containing program options")
		("archive-type", boost::program_options::value<unsigned>(&OPT::nArchiveType)->default_value(2), "project archive type: 0-text, 1-binary, 2-compressed binary")
		("process-priority", boost::program_options::value<int>(&OPT::nProcessPriority)->default_value(-1), "process priority (below normal by default)")
		("max-threads", boost::program_options::value<unsigned>(&OPT::nMaxThreads)->default_value(0), "maximum number of threads (0 for using all available cores)")
		#if TD_VERBOSE != TD_VERBOSE_OFF
		("verbosity,v", boost::program_options::value<int>(&g_nVerbosityLevel)->default_value(
			#if TD_VERBOSE == TD_VERBOSE_DEBUG
			3
			#else
			2
			#endif
			), "verbosity level")
		#endif
		;

	// group of options allowed both on command line and in config file
	boost::program_options::options_description config("Refine options");
	config.add_options()
		("input-file,i", boost::program_options::value<std::string>(&OPT::strInputFileName), "input filename containing camera poses and image list")
		("output-file,o", boost::program_options::value<std::string>(&OPT::strOutputFileName), "output filename for storing the mesh")
		("resolution-level", boost::program_options::value<unsigned>(&OPT::nResolutionLevel)->default_value(0), "how many times to scale down the images before mesh refinement")
		("min-resolution", boost::program_options::value<unsigned>(&OPT::nMinResolution)->default_value(640), "do not scale images lower than this resolution")
		("max-views", boost::program_options::value<unsigned>(&OPT::nMaxViews)->default_value(8), "maximum number of neighbor images used to refine the mesh")
		("decimate", boost::program_options::value<float>(&OPT::fDecimateMesh)->default_value(0.f), "decimation factor in range [0..1] to be applied to the input surface before refinement (0 - auto, 1 - disabled)")
		("close-holes", boost::program_options::value<unsigned>(&OPT::nCloseHoles)->default_value(30), "try to close small holes in the input surface (0 - disabled)")
		("ensure-edge-size", boost::program_options::value<unsigned>(&OPT::nEnsureEdgeSize)->default_value(1), "ensure edge size and improve vertex valence of the input surface (0 - disabled, 1 - auto, 2 - force)")
		("max-face-area", boost::program_options::value<unsigned>(&OPT::nMaxFaceArea)->default_value(64), "maximum face area projected in any pair of images that is not subdivided (0 - disabled)")
		("scales", boost::program_options::value<unsigned>(&OPT::nScales)->default_value(2), "how many iterations to run mesh optimization on multi-scale images")
		("scale-step", boost::program_options::value<float>(&OPT::fScaleStep)->default_value(0.5f), "image scale factor used at each mesh optimization step")
		("reduce-memory", boost::program_options::value<unsigned>(&OPT::nReduceMemory)->default_value(1), "recompute some data in order to reduce memory requirements")
		("alternate-pair", boost::program_options::value<unsigned>(&OPT::nAlternatePair)->default_value(0), "refine mesh using an image pair alternatively as reference (0 - both, 1 - alternate, 2 - only left, 3 - only right)")
		("regularity-weight", boost::program_options::value<float>(&OPT::fRegularityWeight)->default_value(0.2f), "scalar regularity weight to balance between photo-consistency and regularization terms during mesh optimization")
		("rigidity-elasticity-ratio", boost::program_options::value<float>(&OPT::fRatioRigidityElasticity)->default_value(0.9f), "scalar ratio used to compute the regularity gradient as a combination of rigidity and elasticity")
		("gradient-step", boost::program_options::value<float>(&OPT::fGradientStep)->default_value(45.05f), "gradient step to be used instead (0 - auto)")
		("planar-vertex-ratio", boost::program_options::value<float>(&OPT::fPlanarVertexRatio)->default_value(0.f), "threshold used to remove vertices on planar patches (0 - disabled)")
		#ifdef _USE_CUDA
		("use-cuda", boost::program_options::value<bool>(&OPT::bUseCUDA)->default_value(true), "refine mesh using CUDA")
		#endif
		;

	// hidden options, allowed both on command line and
	// in config file, but will not be shown to the user
	boost::program_options::options_description hidden("Hidden options");
	hidden.add_options()
		("mesh-file", boost::program_options::value<std::string>(&OPT::strMeshFileName), "mesh file name to refine (overwrite the existing mesh)")
		;

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(config).add(hidden);

	boost::program_options::options_description config_file_options;
	config_file_options.add(config).add(hidden);

	boost::program_options::positional_options_description p;
	p.add("input-file", -1);

	try {
		// parse command line options
		boost::program_options::store(boost::program_options::command_line_parser((int)argc, argv).options(cmdline_options).positional(p).run(), OPT::vm);
		boost::program_options::notify(OPT::vm);
		INIT_WORKING_FOLDER;
		// parse configuration file
		std::ifstream ifs(MAKE_PATH_SAFE(OPT::strConfigFileName));
		if (ifs) {
			boost::program_options::store(parse_config_file(ifs, config_file_options), OPT::vm);
			boost::program_options::notify(OPT::vm);
		}
	}
	catch (const std::exception& e) {
		LOG(e.what());
		return false;
	}

	// initialize the log file
	OPEN_LOGFILE(MAKE_PATH(APPNAME _T("-")+Util::getUniqueName(0)+_T(".log")).c_str());

	// print application details: version and command line
	Util::LogBuild();
	LOG(_T("Command line:%s"), Util::CommandLineToString(argc, argv).c_str());

	// validate input
	Util::ensureValidPath(OPT::strInputFileName);
	Util::ensureUnifySlash(OPT::strInputFileName);
	if (OPT::vm.count("help") || OPT::strInputFileName.IsEmpty()) {
		boost::program_options::options_description visible("Available options");
		visible.add(generic).add(config);
		GET_LOG() << visible;
	}
	if (OPT::strInputFileName.IsEmpty())
		return false;

	// initialize optional options
	Util::ensureValidPath(OPT::strOutputFileName);
	Util::ensureUnifySlash(OPT::strOutputFileName);
	if (OPT::strOutputFileName.IsEmpty())
		OPT::strOutputFileName = Util::getFullFileName(OPT::strInputFileName) + _T(".mvs");

	// initialize global options
	Process::setCurrentProcessPriority((Process::Priority)OPT::nProcessPriority);
	#ifdef _USE_OPENMP
	if (OPT::nMaxThreads != 0)
		omp_set_num_threads(OPT::nMaxThreads);
	#endif

	#ifdef _USE_BREAKPAD
	// start memory dumper
	MiniDumper::Create(APPNAME, WORKING_FOLDER);
	#endif
	return true;
}
// initialize and parse the command line parameters
bool Initialize_TextureMesh()
{
	size_t argc=7;
	const char * argv[7];
	//for (int i=0;i<9;i++) argv[i]=new char[100];
	argv[0]=(I23D_BIN+"MVSAllSteps").c_str();
	argv[1]="-i";
	argv[2]=(output_dir+"/scene_dense_mesh_refine.mvs").c_str();
	argv[3]="-w";
	argv[4]=working_dir.c_str();
	argv[5]="--resolution-level";
	argv[6]="1";

	// initialize log and console
	OPEN_LOG();
	OPEN_LOGCONSOLE();

	// group of options allowed only on command line
	boost::program_options::options_description generic("Generic options");
	generic.add_options()
		("help,h", "produce this help message")
		("working-folder,w", boost::program_options::value<std::string>(&WORKING_FOLDER), "working directory (default current directory)")
		("config-file,c", boost::program_options::value<std::string>(&OPT::strConfigFileName)->default_value(APPNAME _T(".cfg")), "file name containing program options")
		("archive-type", boost::program_options::value<unsigned>(&OPT::nArchiveType)->default_value(2), "project archive type: 0-text, 1-binary, 2-compressed binary")
		("process-priority", boost::program_options::value<int>(&OPT::nProcessPriority)->default_value(-1), "process priority (below normal by default)")
		("max-threads", boost::program_options::value<unsigned>(&OPT::nMaxThreads)->default_value(0), "maximum number of threads (0 for using all available cores)")
		#if TD_VERBOSE != TD_VERBOSE_OFF
		("verbosity,v", boost::program_options::value<int>(&g_nVerbosityLevel)->default_value(
			#if TD_VERBOSE == TD_VERBOSE_DEBUG
			3
			#else
			2
			#endif
			), "verbosity level")
		#endif
		;

	// group of options allowed both on command line and in config file
	boost::program_options::options_description config("Refine options");
	config.add_options()
		("input-file,i", boost::program_options::value<std::string>(&OPT::strInputFileName), "input filename containing camera poses and image list")
		("output-file,o", boost::program_options::value<std::string>(&OPT::strOutputFileName), "output filename for storing the mesh")
		("resolution-level", boost::program_options::value<unsigned>(&OPT::nResolutionLevel)->default_value(0), "how many times to scale down the images before mesh refinement")
		("min-resolution", boost::program_options::value<unsigned>(&OPT::nMinResolution)->default_value(640), "do not scale images lower than this resolution")
		("outlier-threshold", boost::program_options::value<float>(&OPT::fOutlierThreshold)->default_value(6e-2f), "threshold used to find and remove outlier face textures (0 - disabled)")
		("cost-smoothness-ratio", boost::program_options::value<float>(&OPT::fRatioDataSmoothness)->default_value(0.1f), "ratio used to adjust the preference for more compact patches (1 - best quality/worst compactness, ~0 - worst quality/best compactness)")
		("global-seam-leveling", boost::program_options::value<bool>(&OPT::bGlobalSeamLeveling)->default_value(true), "generate uniform texture patches using global seam leveling")
		("local-seam-leveling", boost::program_options::value<bool>(&OPT::bLocalSeamLeveling)->default_value(true), "generate uniform texture patch borders using local seam leveling")
		("texture-size-multiple", boost::program_options::value<unsigned>(&OPT::nTextureSizeMultiple)->default_value(0), "texture size should be a multiple of this value (0 - power of two)")
		("patch-packing-heuristic", boost::program_options::value<unsigned>(&OPT::nRectPackingHeuristic)->default_value(3), "specify the heuristic used when deciding where to place a new patch (0 - best fit, 3 - good speed, 100 - best speed)")
		("empty-color", boost::program_options::value<uint32_t>(&OPT::nColEmpty)->default_value(0x00FF7F27), "color used for faces not covered by any image")
		;

	// hidden options, allowed both on command line and
	// in config file, but will not be shown to the user
	boost::program_options::options_description hidden("Hidden options");
	hidden.add_options()
		("mesh-file", boost::program_options::value<std::string>(&OPT::strMeshFileName), "mesh file name to texture (overwrite the existing mesh)")
		;

	boost::program_options::options_description cmdline_options;
	cmdline_options.add(generic).add(config).add(hidden);

	boost::program_options::options_description config_file_options;
	config_file_options.add(config).add(hidden);

	boost::program_options::positional_options_description p;
	p.add("input-file", -1);

	try {
		// parse command line options
		boost::program_options::store(boost::program_options::command_line_parser((int)argc, argv).options(cmdline_options).positional(p).run(), OPT::vm);
		boost::program_options::notify(OPT::vm);
		INIT_WORKING_FOLDER;
		// parse configuration file
		std::ifstream ifs(MAKE_PATH_SAFE(OPT::strConfigFileName));
		if (ifs) {
			boost::program_options::store(parse_config_file(ifs, config_file_options), OPT::vm);
			boost::program_options::notify(OPT::vm);
		}
	}
	catch (const std::exception& e) {
		LOG(e.what());
		return false;
	}

	// initialize the log file
	OPEN_LOGFILE(MAKE_PATH(APPNAME _T("-")+Util::getUniqueName(0)+_T(".log")).c_str());

	// print application details: version and command line
	Util::LogBuild();
	LOG(_T("Command line:%s"), Util::CommandLineToString(argc, argv).c_str());

	// validate input
	Util::ensureValidPath(OPT::strInputFileName);
	Util::ensureUnifySlash(OPT::strInputFileName);
	if (OPT::vm.count("help") || OPT::strInputFileName.IsEmpty()) {
		boost::program_options::options_description visible("Available options");
		visible.add(generic).add(config);
		GET_LOG() << visible;
	}
	if (OPT::strInputFileName.IsEmpty())
		return false;

	// initialize optional options
	Util::ensureValidPath(OPT::strOutputFileName);
	Util::ensureUnifySlash(OPT::strOutputFileName);
	if (OPT::strOutputFileName.IsEmpty())
		OPT::strOutputFileName = Util::getFullFileName(OPT::strInputFileName) + _T(".mvs");

	// initialize global options
	Process::setCurrentProcessPriority((Process::Priority)OPT::nProcessPriority);
	#ifdef _USE_OPENMP
	if (OPT::nMaxThreads != 0)
		omp_set_num_threads(OPT::nMaxThreads);
	#endif

	#ifdef _USE_BREAKPAD
	// start memory dumper
	MiniDumper::Create(APPNAME, WORKING_FOLDER);
	#endif
	return true;
}

// finalize application instance
void Finalize()
{
	#if TD_VERBOSE != TD_VERBOSE_OFF
	// print memory statistics
	Util::LogMemoryInfo();
	#endif

	CLOSE_LOGFILE();
	CLOSE_LOGCONSOLE();
	CLOSE_LOG();
}
int MVSRun(std::string Input_dir,std::string Output_dir)
{
	//set< uint32_t > imagesID;
	//for (auto it:images_idx) imagesID.insert(it);
	sfm_data_path=Input_dir+"/sfm_data.json";//std::string(argv[0])+"sfm_data.json";
	output_dir=Output_dir;//=std::string(argv[1]);
	working_dir=output_dir+"/intermediate";
	//MVSAllSteps input_dir output_dir
	#ifdef _DEBUGINFO
	// set _crtBreakAlloc index to stop in <dbgheap.c> at allocation
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);// | _CRTDBG_CHECK_ALWAYS_DF);
	#endif

	if (!Initialize_InterfaceI23dSFM())
		return EXIT_FAILURE;

	MVS::Scene scene(OPT::nMaxThreads);

	TD_TIMER_START();

	if (OPT::bI23D2I23dSFM) {
		// read I23D input data
		if (!scene.Load(MAKE_PATH_SAFE(OPT::strInputFileName)))
			return EXIT_FAILURE;

		// convert data from I23D to I23dSFM
		i23D::MVS_IO::SfM_Scene sceneBAF;
		FOREACH(p, scene.platforms) {
			const MVS::Platform& platform = scene.platforms[p];
			if (platform.cameras.GetSize() != 1) {
				LOG("error: unsupported scene structure");
				return EXIT_FAILURE;
			}
			const MVS::Platform::Camera& camera = platform.cameras[0];
			i23D::MVS_IO::Camera cameraBAF;
			cameraBAF.K = camera.K;
			sceneBAF.cameras.push_back(cameraBAF);
		}
		FOREACH(i, scene.images) {
			const MVS::Image& image = scene.images[i];
			const MVS::Platform& platform = scene.platforms[image.platformID];
			const MVS::Platform::Pose& pose = platform.poses[image.poseID];
			i23D::MVS_IO::Image imageBAF;
			imageBAF.name = image.name;
			imageBAF.name = MAKE_PATH_REL(WORKING_FOLDER_FULL, imageBAF.name);
			imageBAF.id_camera = image.platformID;
			imageBAF.id_pose = (uint32_t)sceneBAF.poses.size();
			sceneBAF.images.push_back(imageBAF);
			i23D::MVS_IO::Pose poseBAF;
			poseBAF.R = pose.R;
			poseBAF.C = pose.C;
			sceneBAF.poses.push_back(poseBAF);
		}
		sceneBAF.vertices.reserve(scene.pointcloud.points.GetSize());
		FOREACH(p, scene.pointcloud.points) {
			const MVS::PointCloud::Point& point = scene.pointcloud.points[p];
			i23D::MVS_IO::Vertex vertexBAF;
			vertexBAF.X = ((const MVS::PointCloud::Point::EVec)point).cast<REAL>();
			const MVS::PointCloud::ViewArr& views = scene.pointcloud.pointViews[p];
			FOREACH(v, views) {
				unsigned viewBAF = views[(uint32_t)v];
				vertexBAF.views.push_back(viewBAF);
			}
			sceneBAF.vertices.push_back(vertexBAF);
		}

		// write I23dSFM input data
		const String strOutputFileNameMVG(OPT::strOutputFileName + MVG_EXT);
		i23D::MVS_IO::ExportScene(MAKE_PATH_SAFE(OPT::strListFileName), MAKE_PATH_SAFE(strOutputFileNameMVG), sceneBAF);

		VERBOSE("Input data exported: %u cameras & %u poses & %u images & %u vertices (%s)", sceneBAF.cameras.size(), sceneBAF.poses.size(), sceneBAF.images.size(), sceneBAF.vertices.size(), TD_TIMER_GET_FMT().c_str());
	} else {
		// convert data from I23dSFM to I23D

		size_t nCameras(0), nPoses(0);

	#ifdef _USE_I23DSFM
	if (OPT::bI23dSFMjson) {
		// read I23dSFM input data from a JSON file
		using namespace i23dSFM::sfm;
		using namespace i23dSFM::cameras;
		SfM_Data sfm_data;
		const String strSfM_Data_Filename(MAKE_PATH_SAFE(OPT::strInputFileName));
		if (!Load(sfm_data, strSfM_Data_Filename, ESfM_Data(ALL))) {
			VERBOSE("error: the input SfM_Data file '%s' cannot be read", strSfM_Data_Filename.c_str());
			return EXIT_FAILURE;
		}
		VERBOSE("Imported data: %u cameras, %u poses, %u images, %u vertices",
				sfm_data.GetIntrinsics().size(),
				sfm_data.GetPoses().size(),
				sfm_data.GetViews().size(),
				sfm_data.GetLandmarks().size());

		// I23dSFM can have not contiguous index, use a map to create the required I23D contiguous ID index
		std::map<i23dSFM::IndexT, uint32_t> map_intrinsic, map_view;

		// define a platform with all the intrinsic group
		nCameras = sfm_data.GetIntrinsics().size();
		for (const auto& intrinsic: sfm_data.GetIntrinsics()) {
			if (isPinhole(intrinsic.second.get()->getType())) {
				const Pinhole_Intrinsic * cam = dynamic_cast<const Pinhole_Intrinsic*>(intrinsic.second.get());
				if (map_intrinsic.count(intrinsic.first) == 0)
					map_intrinsic.insert(std::make_pair(intrinsic.first, scene.platforms.GetSize()));
				MVS::Platform& platform = scene.platforms.AddEmpty();
				// add the camera
				MVS::Platform::Camera& camera = platform.cameras.AddEmpty();
				camera.K = cam->K();
				// sub-pose
				camera.R = RMatrix::IDENTITY;
				camera.C = CMatrix::ZERO;
			}
		}

		// define images & poses
		Util::Progress progress(_T("Processed images"), sfm_data.GetViews().size());
		scene.images.Reserve((uint32_t)sfm_data.GetViews().size());
		for (const auto& view : sfm_data.GetViews()) {
			map_view[view.first] = scene.images.GetSize();
			MVS::Image& image = scene.images.AddEmpty();
			image.name = OPT::strOutputImageFolder + view.second->s_Img_path;

			Util::ensureUnifySlash(image.name);
			image.name = MAKE_PATH_FULL(WORKING_FOLDER_FULL, image.name);
			Util::ensureDirectory(image.name);
			image.platformID = map_intrinsic.at(view.second->id_intrinsic);
			MVS::Platform& platform = scene.platforms[image.platformID];
			image.cameraID = 0;

			i23dSFM::image::Image<i23dSFM::image::RGBColor> imageRGB, imageRGB_ud;
			String pathRoot(sfm_data.s_root_path); Util::ensureDirectorySlash(pathRoot);
			const String srcImage(MAKE_PATH_FULL(WORKING_FOLDER_FULL, pathRoot+view.second->s_Img_path));
			const String& dstImage(image.name);

			if (sfm_data.IsPoseAndIntrinsicDefined(view.second.get())) {
				image.poseID = platform.poses.GetSize();
				MVS::Platform::Pose& pose = platform.poses.AddEmpty();
				const i23dSFM::geometry::Pose3 poseMVG = sfm_data.GetPoseOrDie(view.second.get());
				pose.R = poseMVG.rotation();
				pose.C = poseMVG.center();
				// export undistorted images
				const i23dSFM::cameras::IntrinsicBase * cam = sfm_data.GetIntrinsics().at(view.second->id_intrinsic).get();
				if (cam->have_disto())  {
					// undistort and save the image
					i23dSFM::image::ReadImage(srcImage, &imageRGB);
					i23dSFM::cameras::UndistortImage(imageRGB, cam, imageRGB_ud, i23dSFM::image::BLACK);
					i23dSFM::image::WriteImage(dstImage, imageRGB_ud);
				} else  {
					// no distortion, copy the image
					File::copyFile(srcImage, dstImage);
				}
				++nPoses;
			} else {
				// image have not valid pose, so set an undefined pose
				image.poseID = NO_ID;
				// just copy the image
				File::copyFile(srcImage, dstImage);
			}
			progress.display(scene.images.GetSize());
		}
		progress.close();

		// define structure
		scene.pointcloud.points.Reserve(sfm_data.GetLandmarks().size());
		scene.pointcloud.pointViews.Reserve(sfm_data.GetLandmarks().size());
		for (const auto& vertex: sfm_data.GetLandmarks()) {
			const Landmark & landmark = vertex.second;
			MVS::PointCloud::Point& point = scene.pointcloud.points.AddEmpty();
			point = landmark.X.cast<float>();
			MVS::PointCloud::ViewArr& views = scene.pointcloud.pointViews.AddEmpty();
			for (const auto& observation: landmark.obs)
				views.InsertSort(map_view.at(observation.first));
		}
		
	} else

	#endif
	{
		// read I23dSFM input data from BAF file
		i23D::MVS_IO::SfM_Scene sceneBAF;
		if (!i23D::MVS_IO::ImportScene(MAKE_PATH_SAFE(OPT::strListFileName), MAKE_PATH_SAFE(OPT::strInputFileName), sceneBAF))
			return EXIT_FAILURE;

		// convert data from I23dSFM to I23D
		nCameras = sceneBAF.cameras.size();
		scene.platforms.Reserve((uint32_t)nCameras);
		for (const auto& cameraBAF: sceneBAF.cameras) {
			MVS::Platform& platform = scene.platforms.AddEmpty();
			MVS::Platform::Camera& camera = platform.cameras.AddEmpty();
			camera.K = cameraBAF.K;
			camera.R = RMatrix::IDENTITY;
			camera.C = CMatrix::ZERO;
		}
		nPoses = sceneBAF.images.size();
		scene.images.Reserve((uint32_t)nPoses);
		for (const auto& imageBAF: sceneBAF.images) {
			i23D::MVS_IO::Pose& poseBAF = sceneBAF.poses[imageBAF.id_pose];
			MVS::Image& image = scene.images.AddEmpty();
			image.name = imageBAF.name;
			Util::ensureUnifySlash(image.name);
			image.name = MAKE_PATH_FULL(WORKING_FOLDER_FULL, image.name);
			image.platformID = imageBAF.id_camera;
			MVS::Platform& platform = scene.platforms[image.platformID];
			image.cameraID = 0;
			image.poseID = platform.poses.GetSize();
			MVS::Platform::Pose& pose = platform.poses.AddEmpty();
			pose.R = poseBAF.R;
			pose.C = poseBAF.C;
		}
		scene.pointcloud.points.Reserve(sceneBAF.vertices.size());
		scene.pointcloud.pointViews.Reserve(sceneBAF.vertices.size());
		for (const auto& vertexBAF: sceneBAF.vertices) {
			MVS::PointCloud::Point& point = scene.pointcloud.points.AddEmpty();
			point = vertexBAF.X.cast<float>();
			MVS::PointCloud::ViewArr& views = scene.pointcloud.pointViews.AddEmpty();
			for (const auto& viewBAF: vertexBAF.views)
				views.InsertSort(viewBAF);
		}
	}

		// read images meta-data
		FOREACHPTR(pImage, scene.images) {
			std::cout<<pImage->name<<std::endl;
			IMAGEPTR p(pImage->ReadImageHeader("/edata/lyb/data3/recon/0/mvs/intermediate/undistorted_images/IMG_1869.JPG"));
			std::cout<<p->GetWidth()<<' '<<p->GetHeight()<<' '<<p->GetFileName()<<std::endl;
			//std::cout<<p->GetDataWidth()<<' '<<p->GetDataHeight()<<' '<<p->GetFileName()<<std::endl;

			if (!pImage->ReloadImage(0, false))
				LOG("error: can not read image %s", pImage->name.c_str());
		}
			
		if (OPT::bNormalizeIntrinsics) {
			// normalize camera intrinsics
			FOREACH(p, scene.platforms) {
				MVS::Platform& platform = scene.platforms[p];
				FOREACH(c, platform.cameras) {
					MVS::Platform::Camera& camera = platform.cameras[c];
					// find one image using this camera
					MVS::Image* pImage(NULL);
					FOREACHPTR(pImg, scene.images) {
						if (pImg->platformID == p && pImg->cameraID == c) {
							pImage = pImg;
							break;
						}
					}
					if (pImage == NULL) {
						LOG("error: no image using camera %u of platform %u", c, p);
						continue;
					}
					const REAL fScale(REAL(1)/MVS::Camera::GetNormalizationScale(pImage->width, pImage->height));
					camera.K(0,0) *= fScale;
					camera.K(1,1) *= fScale;
					camera.K(0,2) *= fScale;
					camera.K(1,2) *= fScale;
				}
			}
		}

		// write I23D input data
		scene.Save(MAKE_PATH_SAFE(OPT::strOutputFileName), (ARCHIVE_TYPE)OPT::nArchiveType);

		VERBOSE("Exported data: %u platforms, %u cameras, %u poses, %u images, %u vertices (%s)",
				scene.platforms.GetSize(), nCameras, nPoses, scene.images.GetSize(), scene.pointcloud.GetSize(),
				TD_TIMER_GET_FMT().c_str());
	}
	//DensifyPointCloud
	if (!Initialize_DensifyPointCloud())
		return EXIT_FAILURE;
	//std::cout<<"---------    "<<OPT::strInputFileName<<std::endl;
	if (!scene.Load(MAKE_PATH_SAFE(OPT::strInputFileName)))
		return EXIT_FAILURE;

	FOREACH(idxImage, scene.images) {
		Image& imageData = scene.images[idxImage];
		std::cout<<imageData.width<<' '<<imageData.height<<' '<<imageData.platformID<<' '<<imageData.cameraID<<' '<<imageData.poseID<<' '<<imageData.name<<std::endl;
	}
	
	if (scene.pointcloud.IsEmpty()) {
		VERBOSE("error: empty initial point-cloud");
		return EXIT_FAILURE;
	}
	//TD_TIMER_START();
	if (!scene.DenseReconstruction())
		return EXIT_FAILURE;
	VERBOSE("Densifying point-cloud completed: %u points (%s)", scene.pointcloud.points.GetSize(), TD_TIMER_GET_FMT().c_str());

	String baseFileName(MAKE_PATH_SAFE(Util::getFullFileName(OPT::strOutputFileName) + _T("_dense")));
	scene.Save(baseFileName+_T(".mvs"), (ARCHIVE_TYPE)OPT::nArchiveType);
	scene.pointcloud.Save(baseFileName+_T(".ply"));

	//reconstructMesh===============================================================================================================================
	if (!Initialize_ReconstructMesh())
		return EXIT_FAILURE;
	if (OPT::bMeshExport) {
		// load project
		if (!scene.Load(MAKE_PATH_SAFE(OPT::strInputFileName)) || scene.mesh.IsEmpty())
			return EXIT_FAILURE;
		// save mesh
		const String fileName(MAKE_PATH_SAFE(OPT::strOutputFileName));
		scene.mesh.Save(fileName);
		#if TD_VERBOSE != TD_VERBOSE_OFF
		if (VERBOSITY_LEVEL > 2)
			scene.ExportCamerasMLP(Util::getFullFileName(OPT::strOutputFileName)+_T(".mlp"), fileName);
		#endif
	} else {
		if (OPT::strMeshFileName.IsEmpty()) {
			// load point-cloud and reconstruct a coarse mesh
			if (!scene.Load(MAKE_PATH_SAFE(OPT::strInputFileName)))
				return EXIT_FAILURE;
			// make sure the image neighbors are initialized before deleting the point-cloud
			FOREACH(idxImage, scene.images) {
				const Image& imageData = scene.images[idxImage];
				if (!imageData.IsValid())
					continue;
				if (imageData.neighbors.IsEmpty()) {
					IndexArr points;
					scene.SelectNeighborViews(idxImage, points);
				}
			}
			//TD_TIMER_START();
			if (!scene.ReconstructMesh(OPT::fDistInsert, OPT::bUseFreeSpaceSupport))
				return EXIT_FAILURE;
			VERBOSE("Mesh reconstruction completed: %u vertices, %u faces (%s)", scene.mesh.vertices.GetSize(), scene.mesh.faces.GetSize(), TD_TIMER_GET_FMT().c_str());
			#if TD_VERBOSE != TD_VERBOSE_OFF
			if (VERBOSITY_LEVEL > 2) {
				// dump raw mesh
				scene.mesh.Save(MAKE_PATH_SAFE(Util::getFullFileName(OPT::strOutputFileName) + _T("_mesh_raw.ply")));
			}
			#endif
		} else {
			// load existing mesh to clean
			scene.mesh.Load(MAKE_PATH_SAFE(OPT::strMeshFileName));
		}

		// clean the mesh
		scene.mesh.Clean(OPT::fDecimateMesh, OPT::fRemoveSpurious, OPT::bRemoveSpikes, OPT::nCloseHoles, OPT::nSmoothMesh, false);
		scene.mesh.Clean(1.f, 0.f, OPT::bRemoveSpikes, OPT::nCloseHoles, 0, false); // extra cleaning trying to close more holes
		scene.mesh.Clean(1.f, 0.f, false, 0, 0, true); // extra cleaning to remove non-manifold problems created by closing holes

		// save the final mesh
		baseFileName=(MAKE_PATH_SAFE(Util::getFullFileName(OPT::strOutputFileName) + _T("_mesh")));
		scene.Save(baseFileName+_T(".mvs"), (ARCHIVE_TYPE)OPT::nArchiveType);
		scene.mesh.Save(baseFileName+_T(".ply"));
		#if TD_VERBOSE != TD_VERBOSE_OFF
		if (VERBOSITY_LEVEL > 2)
			scene.ExportCamerasMLP(baseFileName+_T(".mlp"), baseFileName+_T(".ply"));
		#endif
	}

	//refineMesh==============================================================================================================================================================
	if (!Initialize_RefineMesh())
		return EXIT_FAILURE;
	// load and refine the coarse mesh
	if (!scene.Load(MAKE_PATH_SAFE(OPT::strInputFileName)))
		return EXIT_FAILURE;
	if (!OPT::strMeshFileName.IsEmpty()) {
		// load given coarse mesh
		scene.mesh.Load(MAKE_PATH_SAFE(OPT::strMeshFileName));
	}
	if (scene.mesh.IsEmpty()) {
		VERBOSE("error: empty initial mesh");
		return EXIT_FAILURE;
	}
	//TD_TIMER_START();
	#ifdef _USE_CUDA
	if (!OPT::bUseCUDA ||
		!scene.RefineMeshCUDA(OPT::nResolutionLevel, OPT::nMinResolution, OPT::nMaxViews,
							  OPT::fDecimateMesh, OPT::nCloseHoles, OPT::nEnsureEdgeSize,
							  OPT::nMaxFaceArea,
							  OPT::nScales, OPT::fScaleStep,
							  OPT::nAlternatePair>10 ? OPT::nAlternatePair%10 : 0,
							  OPT::fRegularityWeight,
							  OPT::fRatioRigidityElasticity,
							  OPT::fGradientStep))
	#endif
	if (!scene.RefineMesh(OPT::nResolutionLevel, OPT::nMinResolution, OPT::nMaxViews,
						  OPT::fDecimateMesh, OPT::nCloseHoles, OPT::nEnsureEdgeSize,
						  OPT::nMaxFaceArea,
						  OPT::nScales, OPT::fScaleStep,
						  OPT::nReduceMemory, OPT::nAlternatePair,
						  OPT::fRegularityWeight,
						  OPT::fRatioRigidityElasticity,
						  OPT::fPlanarVertexRatio,
						  OPT::fGradientStep))
		return EXIT_FAILURE;
	VERBOSE("Mesh refinement completed: %u vertices, %u faces (%s)", scene.mesh.vertices.GetSize(), scene.mesh.faces.GetSize(), TD_TIMER_GET_FMT().c_str());

	// save the final mesh
	baseFileName=MAKE_PATH_SAFE(Util::getFullFileName(OPT::strOutputFileName) + _T("_refine"));
	scene.Save(baseFileName+_T(".mvs"), (ARCHIVE_TYPE)OPT::nArchiveType);
	scene.mesh.Save(baseFileName+_T(".ply"));
	#if TD_VERBOSE != TD_VERBOSE_OFF
	if (VERBOSITY_LEVEL > 2)
		scene.ExportCamerasMLP(baseFileName+_T(".mlp"), baseFileName+_T(".ply"));
	#endif

	//textureMesh=======================================================================================================================================================
	if (!Initialize_TextureMesh())
		return EXIT_FAILURE;
	// load and texture the mesh
	if (!scene.Load(MAKE_PATH_SAFE(OPT::strInputFileName)))
		return EXIT_FAILURE;
	if (!OPT::strMeshFileName.IsEmpty()) {
		// load given mesh
		scene.mesh.Load(MAKE_PATH_SAFE(OPT::strMeshFileName));
	}
	if (scene.mesh.IsEmpty()) {
		VERBOSE("error: empty initial mesh");
		return EXIT_FAILURE;
	}
	//TD_TIMER_START();
	if (!scene.TextureMesh(OPT::nResolutionLevel, OPT::nMinResolution, OPT::fOutlierThreshold, OPT::fRatioDataSmoothness, OPT::bGlobalSeamLeveling, OPT::bLocalSeamLeveling, OPT::nTextureSizeMultiple, OPT::nRectPackingHeuristic, Pixel8U(OPT::nColEmpty)))
		return EXIT_FAILURE;
	VERBOSE("Mesh texturing completed: %u vertices, %u faces (%s)", scene.mesh.vertices.GetSize(), scene.mesh.faces.GetSize(), TD_TIMER_GET_FMT().c_str());

	// save the final mesh
	baseFileName=MAKE_PATH_SAFE(Util::getFullFileName(OPT::strOutputFileName) + _T("_texture"));
	scene.Save(baseFileName+_T(".mvs"), (ARCHIVE_TYPE)OPT::nArchiveType);
	scene.mesh.Save(baseFileName+_T(".obj"));
	#if TD_VERBOSE != TD_VERBOSE_OFF
	if (VERBOSITY_LEVEL > 2)
		scene.ExportCamerasMLP(baseFileName+_T(".mlp"), baseFileName+_T(".obj"));
	#endif



    Finalize();


	return EXIT_SUCCESS;
}
/*----------------------------------------------------------------*/
