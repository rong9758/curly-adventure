
#pragma once

/*
	Class			: CUser_InterFace_Memory
	
	Company			: HyVision System
	Description		: User Memory
	Version			: 4.0.0.1

	Reporting date	: 2014.10.13
	Name			: Sung - Jin UK

	*
	* Do not modify.
	*
*/
#ifdef _UNICODE
#define tstring std::wstring
#else
#define tstring std::string
#endif

#include "Enum_CMBU.h"
#include "ErrorCode.h"
#include "UserFunc.h"
#include "MyDraw.h" 
#include "HVS_CheckSum.h"
#include <vector>
#include <string>
#include<queue>
#include <functional>
#include <numeric>
#include <set>
#include <algorithm>
#include<unordered_map>

#define	PLATFORM_SITE_COUNT 4
#define IMAGE_CONVERT_INDEX 4
#define IMAGE_CRC_INDEX 8
#define CAPTURE_COUNT	7	
#define ASP_PROCESS		2

/////////////////////////////////
#define INIT_FAIL        1
#define SENSOR_TEMP_FAIL 2
#define MOVECYLINER_FAIL 3
#define AEORCAPTURE_FAIL 4
#define STAGINGTIME_FAIL 5

#define SOFTVERSION    ("1.8.6.0")
#define OFFSETPATH    ("D:\\CurrentOffset\\Tester_MachineSpec.ini")

const int WritePass = 0;
const int OpenRegFail = 1;
const int CreateKeyFail = 2;
const int WriteRegValueFail = 3;

//extern "C" { unsigned int _stdcall WorkTask(void* para); }//import thread pool

//20180424 For Golden
typedef struct _ST_GOLDEN_ITEM
{
	tstring strItemName;
	vector<double> vtGoldenData;
	double dAvgValue;
	void Init()
	{
		dAvgValue = 0.0;
	};
}ST_GoldenItem;

//20180524 For New Golden
typedef struct ST_GOLDENDATA
{
	double dAvgValue;
	double dT0AvgValue;
	tstring strItemName;
	tstring strFailItem;
	tstring strInvalidSN;
	/*set<int> endBarcodeIndex;*/
	vector<int> vtGoldenPara;
	vector<double> vtGoldenData;
	vector<tstring> vtGoldenBarcode;

	ST_GOLDENDATA ()
	{
		dAvgValue = 0.0;
		dT0AvgValue = 0.0;
		strItemName = _T("");
		strFailItem = _T("");
		strInvalidSN = _T("");
		/*endBarcodeIndex.clear();*/
		vtGoldenPara.clear();
		vtGoldenData.clear();
		vtGoldenBarcode.clear();
	};
}ST_GoldenData;

//20180629 For RA CheckData
typedef struct _ST_RA_REVIEW
{
	int iIndex;
	double dData;
	tstring strItemReview;
	void Init()
	{
		iIndex = 0;
		dData = 0.0;
		strItemReview = "";
	}
}ST_RA_REVIEW;

//
template <class T> class find_item_name: public binary_function<T, T, bool>
{
public:
	bool operator()(const T& goldenItem ,const T& cuerrentItem) const
	{
		return goldenItem.strItemName == cuerrentItem.strItemName;
	}
};
template <class T> class calc_avg
{
public:
	void operator()(T& st)
	{
		st.dAvgValue = accumulate(st.vtGoldenData.begin(), st.vtGoldenData.end(), 0.0) / st.vtGoldenData.size();
	}
};
// View&Spec
typedef struct TrendInfo
{
	std::vector<string> m_vtView;
	std::map<string, double>  m_mapLSL;
	std::map<string, double>  m_mapUSL;
	std::map<string, string>  m_mapUnit;

	~TrendInfo()
	{
		m_vtView.clear();
		m_vtView.shrink_to_fit();
		m_mapLSL.clear();
		m_mapUSL.clear();
		m_mapUnit.clear();
	}
	int  GetViewCount(void)
	{
		return m_vtView.size();
	}

	BOOL IsLimited(int i)
	{
		if (m_mapUSL.find(m_vtView[i]) != m_mapUSL.end() 
			|| m_mapLSL.find(m_vtView[i]) != m_mapLSL.end())
		{
			return TRUE;
		}
		return FALSE;
	}

	double GetSpecMin(int i)
	{
		if (m_mapLSL.find(m_vtView[i]) != m_mapLSL.end())
		{
			return m_mapLSL[m_vtView[i]];
		}

		return -65535.0;
	}
	double GetSpecMax(int i)
	{
		if (m_mapUSL.find(m_vtView[i]) != m_mapUSL.end())
		{
			return m_mapUSL[m_vtView[i]];
		}

		return 65535.0;
	}
	string GetUnit(int i)
	{
		if (m_mapUnit.find(m_vtView[i]) != m_mapUnit.end())
		{
			return m_mapUnit[m_vtView[i]];
		}

		return "";
	}
	string GetView(int i)
	{
		if (i < m_vtView.size())
		{
			return m_vtView[i];
		}
		return "N/A";
	}
}TRENDINFO, *PTRENDINFO;
// Algorithm Parameter
typedef struct _VIDEOTYPE 
{
	int     Width;
	int     Height;
	int     nMachine;        //0:Pentagon  1:I300
	PUSHORT *pRawImgSet;     //2 or 5 image set
	char    *szDllPath;      //The path of dll file
	PUCHAR  pOTPRead;        //NVM read map
	PUCHAR  pOTPWrite;       //NVM write map
	PUCHAR  pAspOtpRead;
	int     nNVMStatus;      //0:Empty  1:Not empty
	char    szBarcode[18];   //Barcode from I300
	char    szErrMsg[256];   //Show Error Message to I300
	PTRENDINFO pTrendInfo;
}VIDEOTYPE;

typedef struct MTCPData 
{
	std::vector<std::string> m_vtItem;
	std::vector<std::string> m_vtStarttime;
	std::vector<std::string> m_vtEndTme;
	std::vector<std::string> m_vtTestTime;
	std::vector<std::string> m_vtLow;
	std::vector<std::string> m_vtHight;
	std::vector<std::string> m_vtUnit;
	std::vector<std::string> m_vtValue;

	void InitMTCPData()
	{
		m_vtItem.clear();
		m_vtStarttime.clear();
		m_vtEndTme.clear();
		m_vtTestTime.clear();
		m_vtLow.clear();
		m_vtHight.clear();
		m_vtUnit.clear();
		m_vtValue.clear();
	}

	~MTCPData()
	{
		m_vtItem.clear();
		m_vtStarttime.clear();
		m_vtEndTme.clear();
		m_vtTestTime.clear();
		m_vtLow.clear();
		m_vtHight.clear();
		m_vtUnit.clear();
		m_vtValue.clear();
	}
}MTCPDATA, *PMTCPDATA;

#ifdef DLL_EXPORTS
#define DLL_API __declspec(dllexport)
#else
#define DLL_API __declspec(dllimport)
#endif

struct Res
{
	char* *response_data;
	size_t response_size;

	Res()
	{

	}

	Res(char* *res_data, size_t res_size)
	{
		//memcpy(response_data, res_data, res_size);
		response_data = res_data;
		response_size = res_size;
	}
	//Res(char* *res_data, size_t res_size)
	//{
	//	memcpy(response_data, res_data, res_size);
	//	//response_data = res_data;
	//	response_size = res_size;
	//}
};

struct Config
{
	int maxConfigCount;
	int delayMSeconds;
	int deadlineSeconds;
	int conSeconds;
	int getCount;
};

struct IgRPC_API
{
	//virtual bool Send(vector<const char*> &send_data, vector<CChar> &response_data) = 0;
	//virtual bool Send(vector<const char*> &send_data) = 0;
	virtual bool GetServerIP() = 0;
	//virtual bool GetServerIP_() = 0;
	virtual void SetConfig(Config& conf) = 0;
	virtual void SetLabel(const char* label_c) = 0;
	virtual bool Send(int pivot_type, const char* *send_data, size_t send_data_size) = 0;
	virtual void Release() = 0;
	virtual void Dispose() = 0;
	//virtual Res Response() = 0;
	//virtual void ResponseTest() = 0;
	virtual void Response(Res *res) = 0;
	//virtual void Response(char* *response_data, size_t* response_size) = 0;
	//virtual size_t Response_Size() = 0;
};


extern "C" DLL_API IgRPC_API* __stdcall New();

///220829 remove
////typedef struct TMTCPData
////{
////	std::string info;
////	std::string errCode;
////	std::vector<std::string> m_Group;
////	std::vector<std::string> m_Fom;
////	std::vector<std::string> m_Type;
////	std::vector<std::string> m_Value;
////	std::vector<std::string> m_Note;
////
////	void InitTMTCPDATA()
////	{
////		info = "";
////		errCode = "";
////		m_Group.clear();
////		m_Fom.clear();
////		m_Type.clear();
////		m_Value.clear();
////		m_Note.clear();
////	}
////
////	~TMTCPData()
////	{
////		info = "";
////		errCode = "";
////		m_Group.clear();
////		m_Fom.clear();
////		m_Type.clear();
////		m_Value.clear();
////		m_Note.clear();
////	}
////}TMTCPDATA, *PTMTCPDATA;

//const char MTCPGroup[][32] = {"TSCR", "SIFT", "EE0A", "EE0C", "TMPA", "TMPB", "TMPC", "AE00", "LSCA", "LSCN", "NOIL", "OC00", "RI00", "RU00", "SGNL", "ASPG", "XTLG", "ASPR", "XTLR", "DPCA", "DPCB", "LCB0", "NOID", "LIND", "DSNU", "FPN0", "LINL", "AL60", "SF60", "AL20", "SF20", "PRNU", "STAG", "SHA0", "POST", "TSED"};

//import thread pool - 20220921
//typedef UINT(*Task)(void* lpParam, VIDEOTYPE &, PUSHORT, char *, char *);//import thread pool
typedef UINT(*Task)(void* lpParam, char *);//import thread pool

typedef struct
{
	VIDEOTYPE videotype;
	PUSHORT src;
	char* config_ini;
	char item[128];
	void* lpParam;
}Arg;

class Mutex
{
public:
	Mutex() :_mutex(NULL)
	{
		this->_mutex = CreateMutex(NULL, false, NULL);
	};
	virtual ~Mutex()
	{
		if (this->_mutex)
		{
			CloseHandle(this->_mutex);
		}
		this->_mutex = NULL;
	};

	Mutex(const Mutex& obj) { *this = obj; }

	Mutex& operator=(const Mutex& obj)
	{
		if (this != &obj)
		{
			this->_mutex = obj._mutex;
		}
		return *this;
	}

private:
	HANDLE _mutex;

public:
	void lock()
	{
		WaitForSingleObject(this->_mutex, INFINITE);
	}
	void unlock()
	{
		ReleaseMutex(this->_mutex);
	}

};


class ThreadPool
{

public:
	ThreadPool()
	{
		this->mnMaxNum = 6;
		this->mIsStop = false;
		this->mThreadHandle.resize(this->mnMaxNum);
	};
	~ThreadPool()
	{
		for (int i = 0; i < this->mnMaxNum; i++)
		{
			if (this->mThreadHandle[i] != NULL)
				CloseHandle(this->mThreadHandle[i]);
			this->mThreadHandle[i] = NULL;
		}
	};
public:
	void addTask(const Task task, void* lpParam, char* item)
	{
		this->mMutex.lock();
		this->mTaskQueue.emplace(task);
		//mArg.videotype = videotype; //maybe use
		//mArg.src = src;
		//mArg.config_ini = iniPath;
		strcpy(mArg.item, item);
		mArg.lpParam = lpParam;
		this->mTaskArgQueue.emplace(mArg);
		if (item != NULL)
		{
			RunningItem.insert(std::pair<std::string, bool>(item, true));
		}
			
		this->mMutex.unlock();
	};

	void setMaxThreadNum(int mMaxNum)
	{
		if (mMaxNum	> 0)
		{
			this->mnMaxNum = mMaxNum;
			this->mThreadHandle.resize(this->mnMaxNum);
		}
	};
	void ExitThreadPool()
	{
		this->mIsStop = true;
		for (int i = 0; i < mThreadHandle.size(); i++)
		{
			WaitForSingleObject(mThreadHandle[i], INFINITE);
		}
	};

static unsigned int _stdcall WorkTask(void* para)
	{
		ThreadPool* pool = (ThreadPool*)para;
		while (!pool->mIsStop)
		{
			pool->mMutex.lock();
			if (pool->mTaskQueue.size() <= 0)
			{
				pool->mMutex.unlock();

				Sleep(20);
				continue;
			}
			Task task = pool->mTaskQueue.front();
			Arg nArg = pool->mTaskArgQueue.front();
			pool->mTaskQueue.pop();
			pool->mTaskArgQueue.pop();
			pool->mMutex.unlock();
			task(nArg.lpParam, nArg.item);

			pool->RuningTaskMut.lock();
			std::map<std::string, bool>::iterator it;
			it = pool->RunningItem.find(nArg.item);
			if (it != pool->RunningItem.end())
				pool->RunningItem.erase(it);

			pool->RuningTaskMut.unlock();
		}
		return 0;
	}

	void InitialThreadPool()
	{
		for (int i = 0; i < this->mnMaxNum; i++)
		{
			this->mThreadHandle.at(i) = (HANDLE)_beginthreadex(NULL, 0, WorkTask, this, 0, NULL);
		}
	};

	bool WaitAll(int timeout)
	{
		int count = timeout/10;
		while (count--)
		{
			if (RunningItem.size() <= 0)
			{
				return true;
			}

			Sleep(10);
		}
		if (count <= 0)
			return false;
	};

	bool WaitNoAll(char* str,int timeout)
	{
		int count = timeout/10;
		while (count--)
		{
			WaitItem = RunningItem;
			if (WaitItem.find(str) != WaitItem.end())
			{
				std::map<std::string, bool>::iterator it;
				it = WaitItem.find(str);
				WaitItem.erase(it);
			}
			if (WaitItem.size() <= 0)
			{
				return true;
			}
			Sleep(10);
		}
		if (count <= 0)
			return false;
	};
public:
	Mutex mMutex;
	Mutex RuningTaskMut;
	bool mIsStop;

	queue<Task> mTaskQueue;
	queue<Arg> mTaskArgQueue;
	std::map < std::string, bool> RunningItem;
	std::map<std::string, bool>  WaitItem;
private:
	int mnMaxNum;
	std::vector<HANDLE> mThreadHandle;
	Arg mArg;
};


class CUser_InterFace_Memory
{
public:
	CUserFunc*		pUserFunc;
	CMyDraw*		m_pDrawObj;

	BYTE			*m_pCaptureBuffer[PLATFORM_SITE_COUNT];
	BYTE			*m_pBMPBuffer[PLATFORM_SITE_COUNT];
	BYTE			*m_pYBMPBuffer[PLATFORM_SITE_COUNT];
	BYTE			*m_pRawBuffer[PLATFORM_SITE_COUNT][CAPTURE_COUNT];	//0-4: raw buffer
																		//5 : 2Avg buffer
																		//6 : 5Avg buffer
	WORD			*m_pCropBuffer[PLATFORM_SITE_COUNT][CAPTURE_COUNT];
	WORD			*m_pAPCBuffer[PLATFORM_SITE_COUNT][CAPTURE_COUNT];	//0-4: APC Image Buffer 0 - 4 
	WORD			*m_pWordBuffer[PLATFORM_SITE_COUNT];
										
	
	int				m_nWidth;
	int				m_nHeight;
	VIDEOTYPE		g_VideoType[PLATFORM_SITE_COUNT];

	
	//20180424 For Golden
	/*vector<ST_GoldenItem> m_vtGoldenItem[50];*/
	vector< vector<ST_GoldenItem> > m_vtGoldenItem;

	vector <ST_GoldenData> m_vtGoldenValue;
	UINT m_iGoldenResult;


	CHVSPlatformMemory* m_pMemory;

public:

	int init_flag; // 220707
	//221210
	BOOL Nvm_OTP_Flag ;
	BOOL MTCP_Result_Status ;
	UINT MTCP_ERR_CODE ;
	BOOL Send_MTCP_Flag;

	LPBYTE m_pImageBuffer[12];

//	char m_cSFRcapturePath[65536];

	UINT m_nImageSize;
	UINT m_nImageSize_CRC;
	char m_cArrTestStationNum[4][49];//161002 HVS add
	char sSensorID[4][17]; //161002 HVS add
	char szLensSN[4][20]; // add 20190223
	char szFlexID[4][20]; // add 20190427
	int  m_nInitFlag[4][4];//CMBU add, 20161008, check if init success, [4]:Para, [4]:Site
	int  m_nStopFlag[4];//CMBU add, 20161008, check if Para Stop SiteD, [4]:Para

	int  m_nUploadID[4];//[4]:Para
	int  m_nSiteRet[4];//[4]:Para

	int  m_nParaResult[4];//[4]:Para    2017/08/05
	int  m_bCheckLensFlag[4];//[4]:Para    2017/11/01
	
	BOOL  m_bCheckFlexIDFail[4];//[4]:Para    2add 20190427
	BOOL m_bIsInitFial[4];//[4]:Para 20181207

	BOOL m_bApcFlag[4];//[4]:para
	int  m_iCurRegVal[4]; // [4]:para
	int  m_nNvmStatus[4]; //[4]:Para
	char m_szDefect[4][128];//[4]:Para
// 	char m_szUploadLog[5][102400];  //0:SiteA  1:SiteB  2:SiteC  3:OS  4:RMT
// 
// 	UCHAR m_szTestLog[4][1024]; //[4]:Para
	char *m_szUploadLog[5];  //0:SiteA  1:SiteB  2:SiteC  3:OS  4:RMT
	char *m_szGoldenCheckLog[5];  //0:SiteA  1:SiteB  2:SiteC  3:OS  4:RMT
	char *m_szRawImagePath[5];  //0:SiteA  1:SiteB  2:SiteC  3:OS  4:RMT
	char *m_szUploadHeaderLog[5];  //0:SiteA  1:SiteB  2:SiteC  3:OS  4:RMT
	char *m_szRACheckLog[5];  //0:SiteA  1:SiteB  2:SiteC  3:OS  4:RMT

	char *m_szGoldenAvgData[5];  //0:SiteA  1:SiteB  2:SiteC  3:OS  4:RMT

	UCHAR *m_szTestLog[4]; //[4]:Para

	BOOL  m_bRmtinit;

	PUCHAR m_cBC_ConfigData[4]; 
	double m_dTemperature[4];	
	BOOL   m_bHavePass[4];
	double  m_dsfr60Delta[4][2][3][8];   // Outer&Middle n8&n4&n2 values
	double  m_dsfr20Delta[4][2][3][8];   // Outer&Middle n8&n4&n2 values
	UINT m_iItemCount[4];						//Check Item Counts 20171124

	BOOL   m_bDoubleSNFlag[4];
	BOOL   m_bSFCRet[4];
// 	CString m_strStorePath;

	//20220921 ThreadTool
	ThreadPool* m_ThreadPool[4];
	Mutex* mut[4];
	//ThreadPool* m_ThreadPool;


	//20200720
	std::unordered_map<std::string, MTCPDATA> m_mapMTCPData[4];
	//std::map<std::string, MTCPDATA> m_mapMTCPData[4];
	std::map<std::string, std::string> m_mapTestVer[4];
	std::string m_strOvenTime[4];
	std::string m_strOvenFinishDate[4];

	std::map<tstring, std::vector<tstring>> m_mapSFCInfo;
	std::map<tstring, tstring> m_mapSFCInfo_LensSN;
	PUCHAR  m_pNVMBeforeWriteData[4];
	PUCHAR  m_pNVMAfterWriteData[4];

	char m_szBarcode[4][18]; //20201102 add
	char m_szASICID[4][18]; //20201107 add
	int m_iPRE_FAULT[4];//20201204 add
	int m_iPRE_RECORDER[4];//20201204 add
	int m_bUOPFlag[4];//20201225 add
	std::string m_strLogFilePath[4]; //20201130 add

	CUser_InterFace_Memory()
	{

	}	

	~CUser_InterFace_Memory()
	{
		
	}

public:

	CUserFunc* GetUserFunc()
	{
		return pUserFunc;
	}

	CMyDraw* GetDrawObject()
	{
		return m_pDrawObj;
	}

	BYTE* GetCaptureBuffer(UINT nUsbID)
	{
		return m_pCaptureBuffer[nUsbID];
	}

	BYTE* GetBMPBuffer(UINT nUsbID)
	{
		return m_pBMPBuffer[nUsbID];
	}

	BYTE* GetYBMPBuffer(UINT nUsbID)
	{
		return m_pYBMPBuffer[nUsbID];
	}
	
	USHORT* GetAPCBuffer(UINT nUsbID, int nCaptureBuffer)
	{
		return m_pAPCBuffer[nUsbID][nCaptureBuffer];
	}

	USHORT** GetAPCBuffer(UINT nUsbID)
	{
		return m_pAPCBuffer[nUsbID];
	}

	BYTE* GetRawBuffer(UINT nUsbID, UINT nCaptureBuffer)
	{
		return m_pRawBuffer[nUsbID][nCaptureBuffer];
	}

	WORD* GetCropBuffer(UINT nUsbID, UINT nCaptureBuffer)
	{
		return m_pCropBuffer[nUsbID][nCaptureBuffer];
	}

	WORD** GetCropBuffer(UINT nUSbID)
	{
		return m_pCropBuffer[nUSbID];
	}

	WORD* GetWordBuffer(UINT nUsbID)
	{
		return m_pWordBuffer[nUsbID];
	}

	void CleanImageBuffer(UINT nUsbID, int nWidth, int nHeight)
	{
		for (int i = 0; i < CAPTURE_COUNT; ++i)
		{
			memset(m_pRawBuffer[nUsbID][i], 0, sizeof(BYTE)*nWidth*nHeight * 2);
			memset(m_pCropBuffer[nUsbID][i], 0, sizeof(WORD)*nWidth*nHeight * 2);
			memset(m_pAPCBuffer[nUsbID][i], 0, sizeof(WORD)*nWidth*nHeight * 2);
		}
	}
	//TODO 220829 add test_ver
	bool GetTestVer(int iPara,string test_item,string test_view,string test_value)
	{
		bool test_ver_flag = true;
		//std::transform(test_view.begin(), test_view.end(), test_view.begin(), ::tolower);
		if ("test_ver" == test_view || "Test_ver" == test_view || "sfr_test_ver" == test_view)
		{
			if (m_mapTestVer[iPara].find(test_item) == m_mapTestVer[iPara].end())
			{
				m_mapTestVer[iPara].insert(std::pair<string, string>(test_item, test_value));
			}
			else
			{
				if (m_mapTestVer[iPara].find(test_item) !=m_mapTestVer[iPara].end())
				{
					m_mapTestVer[iPara][test_item] = test_value;
				}
			}
		}
		else
		{
			test_ver_flag = false;
		}
		return test_ver_flag;
	}

	int WriteREG_SZ(LPCTSTR lpKeyFolder, LPCTSTR lpSubKey, LPCTSTR lpValueName, LPCTSTR wcData)
	{
		int iResult = 0;
		HKEY hKey = NULL;			
		HKEY hTempKey = NULL;
		DWORD dwLen = 0;

		dwLen = strlen(wcData) * sizeof(wchar_t);

		if (ERROR_SUCCESS == ::RegOpenKeyEx(HKEY_CURRENT_USER, lpKeyFolder, 0, KEY_SET_VALUE, &hKey))
		{
			if (ERROR_SUCCESS == ::RegCreateKey(hKey, lpSubKey, &hTempKey))
			{
				if (ERROR_SUCCESS == ::RegSetValueEx(hTempKey, lpValueName, 0, REG_SZ, (const BYTE*)wcData, dwLen))
				{
					iResult = WritePass;
				}
				else
				{
					iResult = WriteRegValueFail;
				}
			}
			else
			{
				iResult = CreateKeyFail;
			}
		}
		else
		{
			iResult = OpenRegFail;
		}

		::RegCloseKey(hKey);

		return iResult;
	}

	BOOL ReadREG_SZ(LPCTSTR lpKeyFolder, LPCTSTR lpValueName, TCHAR *wcData)
	{
		BOOL bResult = FALSE;
		HKEY hKey = NULL;		
		TCHAR dwValue[256] = {0};
		DWORD dwType = REG_SZ;
		DWORD dwSize = sizeof(TCHAR)*256;

		if (ERROR_SUCCESS == ::RegOpenKeyEx(HKEY_CURRENT_USER, lpKeyFolder, 0, KEY_READ, &hKey))
		{
			if (ERROR_SUCCESS ==  ::RegQueryValueEx(hKey, lpValueName, 0, &dwType, (LPBYTE)&dwValue, &dwSize))
			{
				wsprintf(wcData, dwValue);
				bResult = TRUE;
			}
			else
			{
				bResult = FALSE;
			}
		}
		else
		{
			bResult = FALSE;
		}

		::RegCloseKey(hKey);

		return bResult;
	}



// 	void CreateListLogDirectory(void)
// 	{
// 		GetStoreListLogPath();
// 		CreateAbsoluteDir(m_strStorePath);
// 	}
// 
// 	void GetStoreListLogPath()
// 	{
// 		char szTemp[256] = {0};
// 		char szConfigPath[256] ={0};
// 		CString strMainPath = "";
// 		char szCaseConfigPath[256] = {0};
// 		SYSTEMTIME sysTime;
// 
// 		char szPath_tmp[256] = {0};
// 		::GetModuleFileName(::AfxGetInstanceHandle(), szPath_tmp, 256);
// 		PathRemoveFileSpec(szPath_tmp);
// 		strMainPath = szPath_tmp;
// 		
// 
// 		wsprintf(szCaseConfigPath, "%s\\Config\\Case.ini", strMainPath);
// 		GetPrivateProfileString("Client", "ConfigPath", NULL, szTemp, sizeof(szTemp), szCaseConfigPath);
// 
// 		/*wsprintf(szConfigPath, "%s\\Config\\config.ini", strMainPath);*/
// 		wsprintf(szConfigPath, "%s\\Case\\%s\\Config.ini", strMainPath, szTemp);
// 
// 		GetPrivateProfileString("Config", "StoreDir", "", szTemp, sizeof(szTemp), szConfigPath);
// 		m_strStorePath = szTemp;
// 		m_strStorePath += "\\";
// 
// 		GetLocalTime(&sysTime);
// 		sprintf(szTemp, "ListLog\\%04d%02d%02d", sysTime.wYear, sysTime.wMonth, sysTime.wDay);
// 		m_strStorePath += szTemp;
// 		m_strStorePath += "\\";
// 	}
// 
// 	void CreateAbsoluteDir(CString path)
// 	{
// 		CString str = "";
// 		char    *pTemp = path.GetBuffer(path.GetLength());
// 		int i = 3;
// 
// 		str += pTemp[0];
// 		str += pTemp[1];
// 		str += pTemp[2];
// 		while (pTemp[i] != '\0')
// 		{
// 			str += pTemp[i];
// 			if (pTemp[i] == '\\')
// 			{
// 				CreateDirectory(str.GetBuffer(str.GetLength()), NULL);
// 				str.ReleaseBuffer();
// 			}
// 			i++;
// 		}
// 		path.ReleaseBuffer();
// 	}
// 
// 	void ListLog(UINT nSite, LPCSTR format,...)
// 	{
// 		CString		sResult = "";
// 
// 		va_list args;
// 		va_start(args,format);
// 		sResult.FormatV(format,args);
// 		va_end(args);
// 
// 		SYSTEMTIME		st;
// 		HANDLE				hFile;
// 		DWORD					dw;
// 		CString					tmpBuff = "";
// 		char						str_file[1024];
// 		CString					path = "";
// 		char						szTestSite[4][64] = {"A", "B", "C", "D"};
// 		GetLocalTime(&st);
// 		path.Format("%s%04d_%02d_%02d_ListLog_%s.txt",m_strStorePath, st.wYear, st.wMonth, st.wDay, szTestSite[nSite]);
// 		wsprintf(str_file, path);
// 		hFile = CreateFile(str_file,GENERIC_WRITE,FILE_SHARE_READ,NULL,
// 			OPEN_ALWAYS/*CREATE_ALWAYS*/,FILE_ATTRIBUTE_NORMAL,NULL);
// 
// 		if( hFile == INVALID_HANDLE_VALUE ){
// 			CloseHandle(hFile);
// 			return ;
// 		}
// 		SetFilePointer(hFile,0,NULL,FILE_END);
// 
// 		tmpBuff.Format("[%02d/%02d %02d:%02d:%02d %03d] %s\r\n", st.wMonth, st.wDay, st.wHour, st.wMinute, st.wSecond, st.wMilliseconds, sResult);
// 		WriteFile(hFile,tmpBuff.GetBuffer(tmpBuff.GetLength()),tmpBuff.GetLength(),&dw,NULL);
// 		tmpBuff.ReleaseBuffer();
// 		CloseHandle(hFile);
// 	}
};

