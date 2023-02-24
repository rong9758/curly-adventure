#include "StdAfx.h"
#include "Test.h"

#include <algorithm>
#include "../../Inc/CMBU/ICheckSum.h"
#ifdef WIN64
#pragma comment(lib, "../../../Inc/lib/CheckSum_64.lib")
#else
#pragma comment(lib, "../../../Inc/lib/CheckSum.lib")
#endif

//CTest m_Test;

CTest::CTest()
	: m_pLoginDlg(NULL)
{
	m_pSpecDlg = NULL;
	m_pUserMemory = NULL;
}

CTest::~CTest()
{
	char szExePath[MAX_PATH] = {0};
	///220823 remove
	////sprintf(szExePath, "%s\\DolanDemo\\KillDemo.bat", m_pPfMemory->GetRoot());
	////ShellExecute(NULL, "open", szExePath, NULL, NULL, SW_HIDE);

	sprintf(szExePath, "%s\\KillMTCPTool.bat", m_pPfMemory->GetRoot());
	//AfxMessageBox(LPCTSTR("KillMTCPTool"));
	////m_pPfMemory->SetLog(3, "KillMTCPTool: %s", szExePath);

	ShellExecute(NULL, "open", szExePath, NULL, NULL, SW_HIDE);

	if (m_pLoginDlg != NULL)
	{
		m_pLoginDlg->DestroyWindow();
		delete m_pLoginDlg;
		m_pLoginDlg = NULL;
	}

	if (m_pSpecDlg != NULL)
	{
		m_pSpecDlg->DestroyWindow();
		delete m_pSpecDlg;
		m_pSpecDlg = NULL;
	}

	if (m_pUserMemory)
	{
		for (int n = 0; n < 5; n++)
		{
			if(m_pUserMemory->m_szUploadLog[n])	
			{
				delete m_pUserMemory->m_szUploadLog[n];
			}

			if(m_pUserMemory->m_szGoldenCheckLog[n])	
			{
				delete m_pUserMemory->m_szGoldenCheckLog[n];
			}

			if(m_pUserMemory->m_szGoldenAvgData[n])	
			{
				delete m_pUserMemory->m_szGoldenAvgData[n];
			}

			if(m_pUserMemory->m_szRawImagePath[n])	
			{
				delete m_pUserMemory->m_szRawImagePath[n];
			}

			if(m_pUserMemory->m_szUploadHeaderLog[n])	
			{
				delete m_pUserMemory->m_szUploadHeaderLog[n];
			}
			if(m_pUserMemory->m_szRACheckLog[n])
			{
				delete m_pUserMemory->m_szRACheckLog[n];
			}
		}

		for(int i = 0; i < 4; i++)
		{
			if(m_pUserMemory->m_pCaptureBuffer[i])	delete m_pUserMemory->m_pCaptureBuffer[i];
			if(m_pUserMemory->m_pBMPBuffer[i])		delete m_pUserMemory->m_pBMPBuffer[i];
			if (m_pUserMemory->m_pYBMPBuffer[i])		delete m_pUserMemory->m_pYBMPBuffer[i];
			if(m_pUserMemory->m_pWordBuffer[i])		delete m_pUserMemory->m_pWordBuffer[i];

			if(m_pUserMemory->m_szTestLog[i])	delete m_pUserMemory->m_szTestLog[i];

			for(int j = 0; j < CAPTURE_COUNT; j++)
			{
				if(m_pUserMemory->m_pRawBuffer[i][j])				 delete m_pUserMemory->m_pRawBuffer[i][j];
				if(m_pUserMemory->m_pCropBuffer[i][j])				 delete m_pUserMemory->m_pCropBuffer[i][j];
				if(m_pUserMemory->m_pAPCBuffer[i][j])				 delete m_pUserMemory->m_pAPCBuffer[i][j];		
			}			
			for (int j = 0; j < 5; j++)
			{
				delete m_pUserMemory->g_VideoType[i].pRawImgSet[j];
			}
			delete m_pUserMemory->g_VideoType[i].pRawImgSet;
			
			delete m_pUserMemory->g_VideoType[i].pOTPRead;
			delete m_pUserMemory->g_VideoType[i].pOTPWrite;
			delete m_pUserMemory->g_VideoType[i].pAspOtpRead;

			if (m_pUserMemory->m_pNVMBeforeWriteData[i])  delete [] m_pUserMemory->m_pNVMBeforeWriteData[i];
			if (m_pUserMemory->m_pNVMAfterWriteData[i])  delete [] m_pUserMemory->m_pNVMAfterWriteData[i];
	
			if (m_pUserMemory->m_ThreadPool[i]) { m_pUserMemory->m_ThreadPool[i]->ExitThreadPool(); delete  m_pUserMemory->m_ThreadPool[i]; }
			if (m_pUserMemory->mut[i]) delete m_pUserMemory->mut[i];
		}

		if (m_pUserMemory->pUserFunc)		delete  m_pUserMemory->pUserFunc;
		if (m_pUserMemory->m_pDrawObj)		delete	m_pUserMemory->m_pDrawObj;

	}
	if(m_pUserMemory)				delete m_pUserMemory;
}

BOOL CTest::OnInit(CHVS_InterFace_SW* hvs_SW, void* hvs_HW, void* hvs_Protocol, CHVSPlatformMemory* pMemory)
{
	char szCaption[256] = {0};
	char szTmpCaption[256] = {0};
	char szTemp[256] = {0};
	char szReleaseTime[256] = {0};
	char szCaseConfigPath[MAX_PATH] ={0};
	HWND hwndPlatform = NULL;

	m_pSoftware = hvs_SW;
	m_pHardware = (CHVS_InterFace_HW_V5FP*)hvs_HW;
	m_pProtocol = hvs_Protocol;
	m_pPfMemory = pMemory;

	if (m_pLoginDlg == NULL)
	{
  		m_pLoginDlg = new CLogin(m_pSoftware, m_pHardware, m_pProtocol, m_pPfMemory);
	}
// 		CLogin aaa(m_pSoftware, m_pHardware, m_pProtocol, m_pPfMemory);
// 		aaa.DoModal();
	if (m_pLoginDlg->DoModal() == IDOK)
	{
	}
	else
	{
		ExitProcess(0);
	}
	// Get Platform main window

	tstring strRetest = "";
// 	BOOL bIsRetest = CheckMachineisRetest();
// 	if (bIsRetest)
// 	{
// 		strRetest = "ReTest";
// 	}
// 	else
// 	{
// 		strRetest = "FirstTest";
// 	}


	hwndPlatform = m_pPfMemory->GetPlatformHwnd();

	wsprintf(szCaseConfigPath, "%s\\Config\\Case.ini", m_pPfMemory->GetRoot());
	GetWindowText(hwndPlatform, szCaption, sizeof(szCaption));

	char* pTail = strchr(szCaption, '-');
	char* pHead = strchr(szCaption, ':');
	*pHead = 0;

	GetPrivateProfileString("SysConfig", "ReleaseTime", "", szReleaseTime, sizeof(szReleaseTime), szCaseConfigPath);
	GetPrivateProfileString("Client", "ConfigPath", "", szTemp, sizeof(szTemp), szCaseConfigPath);
	sprintf(szTmpCaption, "%s: %s)%s_%s", szCaption, szReleaseTime, pTail, szTemp);
	SetWindowText(hwndPlatform, szTmpCaption);

	if (m_pSpecDlg == NULL)
	{
		m_pSpecDlg = new CSpecDlg(m_pSoftware, m_pHardware, m_pProtocol, m_pPfMemory);
		m_pSpecDlg->Create(IDD_DIALOG_SPEC, NULL);
	}

	return 0;
}

UINT CTest::StartTest(UINT nUsbID, UINT nSite, UINT nTestNumber)
{
	SYSTEMTIME st;
	GetLocalTime(&st);

	m_pPfMemory->SetLog(nSite, "Release on 2021/01/09 CheckConfig,LotNo");
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

	for (int iPara = 0; iPara < 4; ++iPara)
	{
		//StartTest Image Delete
		char sfrPath[MAX_PATH] = {0};
		sprintf(sfrPath, "C:\\Temp\\%d", iPara);
		sprintf(sfrPath, "%s\\SFRcaptureimage.raw", sfrPath);
		DeleteFile(sfrPath);
		m_pPfMemory->SetLog(nSite, "Delete : %s", sfrPath);
	}


 	m_pPfMemory->SetLog(nSite, "PowerUpSequence");
 	for (int i = 0; i < 4; ++i)
 	{
 		m_pHardware->PowerUpSequence(i, 0, i, m_pHardware->GetOscType());
 	}
 	
 	Sleep(100);
 
 	m_pPfMemory->SetLog(nSite, "PowerDownSequence");
 	for (int j = 0; j < 4; ++j)
 	{
 		m_pHardware->PowerDownSequence(j, 0);
 	}



	if (m_pUserMemory == NULL)
	{
		m_pUserMemory = new CUser_InterFace_Memory;

		m_pUserMemory->init_flag = 0; // 220707

		//221230 add initial
		m_pUserMemory->Nvm_OTP_Flag = FALSE;
		m_pUserMemory->MTCP_Result_Status = FALSE;
		m_pUserMemory->Send_MTCP_Flag = TRUE;
		m_pUserMemory->MTCP_ERR_CODE = 0;

		m_pUserMemory->m_nImageSize = m_pHardware->GetImageSize(m_pHardware->GetDataFormat(), m_pHardware->GetImageWidth(), m_pHardware->GetImageHeight());
		m_pUserMemory->m_nImageSize_CRC = m_pHardware->GetCRCImageSize(m_pHardware->GetDataFormat(), m_pHardware->GetImageWidth(), m_pHardware->GetImageHeight());

		int nWidth = m_pPfMemory->IniRead_Int("SENSOR_SET", "HVS_SENSOR_SIZE_X", 0, m_pPfMemory->GetSensorPath());
		int nHeight = m_pPfMemory->IniRead_Int("SENSOR_SET", "HVS_SENSOR_SIZE_Y", 0, m_pPfMemory->GetSensorPath());

		for (int n = 0; n < 5; n++)
		{
			m_pUserMemory->m_szUploadLog[n] = new char[102400];
			memset(m_pUserMemory->m_szUploadLog[n], 0, sizeof(char)*102400);

			m_pUserMemory->m_szGoldenCheckLog[n] = new char[10240];
			memset(m_pUserMemory->m_szGoldenCheckLog[n], 0, sizeof(char)*10240);

			m_pUserMemory->m_szGoldenAvgData[n] = new char[10240];
			memset(m_pUserMemory->m_szGoldenAvgData[n], 0, sizeof(char)*10240);

			m_pUserMemory->m_szRawImagePath[n] = new char[20480];
			memset(m_pUserMemory->m_szRawImagePath[n], 0, sizeof(char)*20480);

			m_pUserMemory->m_szUploadHeaderLog[n] = new char[1024*200];
			memset(m_pUserMemory->m_szUploadHeaderLog[n], 0, sizeof(char)*1024*200);

			m_pUserMemory->m_szRACheckLog[n] = new char[10240];
			memset(m_pUserMemory->m_szRACheckLog[n], 0, sizeof(char)*10240);
		}

		memset(m_pUserMemory->m_dsfr60Delta, 0, sizeof(m_pUserMemory->m_dsfr60Delta));
		memset(m_pUserMemory->m_dsfr20Delta, 0, sizeof(m_pUserMemory->m_dsfr20Delta));

		for(int i = 0; i < 4; i++)
		{
			m_pUserMemory->m_bSFCRet[i] = TRUE;
			m_pUserMemory->m_bDoubleSNFlag[i] = FALSE;
			m_pUserMemory->m_iItemCount[i] = 0;
			m_pUserMemory->m_bHavePass[i] = FALSE;
			m_pUserMemory->m_bUOPFlag[i] = 0;
			m_pUserMemory->m_pCaptureBuffer[i]	= new BYTE[nWidth*nHeight*3];
			m_pUserMemory->m_pBMPBuffer[i]		= new BYTE[nWidth*nHeight*3];
			m_pUserMemory->m_pYBMPBuffer[i]		= new BYTE[nWidth*nHeight*3];
			m_pUserMemory->m_pWordBuffer[i]		= new WORD[nWidth*nHeight*3];			
			m_pUserMemory->m_cBC_ConfigData[i]  = new UCHAR[OTP_BufSize];

			m_pUserMemory->m_szTestLog[i]	= new UCHAR[1024];

			for(int j = 0 ; j < CAPTURE_COUNT; j++)
			{
				m_pUserMemory->m_pRawBuffer[i][j] = new BYTE[nWidth*nHeight*2];
				m_pUserMemory->m_pCropBuffer[i][j] = new WORD[nWidth*nHeight*2];
				m_pUserMemory->m_pAPCBuffer[i][j] = new WORD[nWidth*nHeight*2];
			}	
			m_pUserMemory->g_VideoType[i].pRawImgSet = new PUSHORT[5];

			for (int j = 0; j < 5; j++)
			{
				m_pUserMemory->g_VideoType[i].pRawImgSet[j] = new USHORT[nWidth * nHeight * 2];
			}

			m_pUserMemory->g_VideoType[i].pOTPRead = NULL;
			m_pUserMemory->g_VideoType[i].pOTPWrite = NULL;
			m_pUserMemory->g_VideoType[i].pAspOtpRead = NULL;

			m_pUserMemory->g_VideoType[i].pOTPWrite  = new UCHAR[NVMSIZE];
			m_pUserMemory->g_VideoType[i].pOTPRead   = new UCHAR[NVMSIZE];		
			m_pUserMemory->g_VideoType[i].pAspOtpRead = new UCHAR[NVMSIZE];

			m_pUserMemory->m_bApcFlag[i] = 0;
			m_pUserMemory->m_iCurRegVal[i] = 0;
			m_pUserMemory->m_nNvmStatus[i] = 0;
			m_pUserMemory->m_dTemperature[i] = 0.0;

			memset(m_pUserMemory->m_pCaptureBuffer[i], 0, sizeof(BYTE)*nWidth*nHeight * 3);
			memset(m_pUserMemory->m_pBMPBuffer[i], 0, sizeof(BYTE)*nWidth*nHeight * 3);
			memset(m_pUserMemory->m_pYBMPBuffer[i], 0, sizeof(BYTE)*nWidth*nHeight * 3);
			memset(m_pUserMemory->m_pWordBuffer[i], 0, sizeof(WORD)*nWidth*nHeight * 3);

			memset(m_pUserMemory->m_szTestLog[i], 0, sizeof(UCHAR)*1024);

			memset(m_pUserMemory->m_cBC_ConfigData[i], 0, OTP_BufSize * sizeof(UCHAR));

			for (int j = 0; j < CAPTURE_COUNT; j++)
			{
				memset(m_pUserMemory->m_pRawBuffer[i][j], 0, sizeof(BYTE)*nWidth*nHeight * 2);
				memset(m_pUserMemory->m_pCropBuffer[i][j], 0, sizeof(WORD)*nWidth*nHeight * 2);
				memset(m_pUserMemory->m_pAPCBuffer[i][j], 0, sizeof(WORD)*nWidth*nHeight * 2);
			}

			for (int j = 0; j < 5; j++)
			{
				memset(m_pUserMemory->g_VideoType[i].pRawImgSet[j], 0, sizeof(USHORT)*nWidth*nHeight*2);
			}

			memset(m_pUserMemory->g_VideoType[i].pOTPWrite, NULL, NVMSIZE * sizeof(UCHAR));
			memset(m_pUserMemory->g_VideoType[i].pOTPRead, NULL, NVMSIZE * sizeof(UCHAR));
			memset(m_pUserMemory->g_VideoType[i].pAspOtpRead, NULL, NVMSIZE * sizeof(UCHAR));
			m_pUserMemory->g_VideoType[i].nNVMStatus = 0;
			memset(m_pUserMemory->g_VideoType[i].szBarcode, NULL, sizeof(m_pUserMemory->g_VideoType->szBarcode));
			memset(m_pUserMemory->g_VideoType[i].szErrMsg, NULL, sizeof(m_pUserMemory->g_VideoType->szErrMsg));

			m_pUserMemory->m_nUploadID[i] = 0;
			m_pUserMemory->m_nSiteRet[i] = 0;

			//2017/08/05
			m_pUserMemory->m_nParaResult[i] = 0;
			m_pUserMemory->m_bCheckLensFlag[i] = 0;
			//20190427
			m_pUserMemory->m_bCheckFlexIDFail[i] = FALSE;

			//2018/12/07
			m_pUserMemory->m_bIsInitFial[i] = FALSE;

			for (int j = 0; j < 4; ++j)
			{
				m_pUserMemory->m_nInitFlag[i][j] = 0;
			}

			m_pUserMemory->m_mapMTCPData[i].clear();
			memset(m_pUserMemory->sSensorID[i], NULL, sizeof(m_pUserMemory->sSensorID[i]));
			memset(m_pUserMemory->m_szBarcode[i], NULL, sizeof(m_pUserMemory->m_szBarcode[i]));
			memset(m_pUserMemory->m_szASICID[i], NULL, sizeof(m_pUserMemory->m_szASICID[i]));
			memset(m_pUserMemory->szFlexID[i], NULL, sizeof(m_pUserMemory->szFlexID[i]));
			memset(m_pUserMemory->szLensSN[i], NULL, sizeof(m_pUserMemory->szLensSN[i]));
// 			memset(m_pUserMemory->m_iPRE_FAULT[i], NULL, sizeof(m_pUserMemory->m_iPRE_FAULT[i]));
// 			memset(m_pUserMemory->m_iPRE_RECORDER[i], NULL, sizeof(m_pUserMemory->m_iPRE_RECORDER[i]));
			m_pUserMemory->m_iPRE_FAULT[i] = 0;
			m_pUserMemory->m_iPRE_RECORDER[i] = 0;
			m_pUserMemory->m_strOvenTime[i] = "";
			m_pUserMemory->m_strOvenFinishDate[i] = "";
			m_pUserMemory->m_strLogFilePath[i].clear();

			m_pUserMemory->m_pNVMBeforeWriteData[i] = new UCHAR[NVMSIZE];
			ZeroMemory(m_pUserMemory->m_pNVMBeforeWriteData[i], NVMSIZE * sizeof(UCHAR));

			m_pUserMemory->m_pNVMAfterWriteData[i] = new UCHAR[NVMSIZE];
			ZeroMemory(m_pUserMemory->m_pNVMAfterWriteData[i], NVMSIZE * sizeof(UCHAR));
		}

		m_pUserMemory->pUserFunc = new CUserFunc(m_pSoftware, m_pHardware, m_pProtocol, m_pPfMemory);
		m_pUserMemory->m_pDrawObj = new CMyDraw;

		for (int j = 0; j < 4; j++)
		{
			m_pUserMemory->m_ThreadPool[j] = new ThreadPool;
			m_pUserMemory->m_ThreadPool[j]->setMaxThreadNum(12);
			m_pUserMemory->m_ThreadPool[j]->InitialThreadPool();
			m_pUserMemory->mut[j] = new Mutex;
		}

		m_pUserMemory->m_bRmtinit = FALSE;
		m_pSoftware->SetUserMemory(m_pUserMemory, "USER_MEMORY");
	}
	m_pPfMemory->RegeditWrite_Int(SHCU, REG_SYSTEM_KEY, "VersionByPass", 1);

	CreateDirectory("D:\\CurrentOffset", NULL);

	char	szMes[80] = {0, };
	int		nOffset_value = 0;
	char strPath[256] = {0};

 	/*wsprintf(strPath, "%s\\Config\\Tester_MachineSpec.ini", m_pPfMemory->GetRoot());*/
	wsprintf(strPath, OFFSETPATH);

	m_pPfMemory->SetLog(0, "strPath: %s", strPath);

	for(int nPara=0; nPara<4; nPara++)
	{
		for(int iCH=0; iCH<6; iCH++)
		{
			wsprintf(szMes, "Current_Dynamic_Measure_offset_%d_%d", iCH, nPara);		
			nOffset_value = m_pPfMemory->IniRead_Int("CURRENT_TEST", szMes, 0, (LPSTR)(LPCTSTR)strPath);
			m_pPfMemory->SetLog(0, "Dynamic %s: %d",szMes, nOffset_value);

			wsprintf(szMes, "Current_Standby_Measure_offset_%d_%d", iCH, nPara);		
			nOffset_value = m_pPfMemory->IniRead_Int("CURRENT_TEST_STANDBY_500uA", szMes, 0, (LPSTR)(LPCTSTR)strPath);
			m_pPfMemory->SetLog(0, "STANDBY_500uA %s: %d",szMes, nOffset_value);

			wsprintf(szMes, "Current_Standby_Measure_offset_%d_%d", iCH, nPara);		
			nOffset_value = m_pPfMemory->IniRead_Int("CURRENT_TEST_STANDBY_5mA", szMes, 0, (LPSTR)(LPCTSTR)strPath);
			m_pPfMemory->SetLog(0, "STANDBY_5mA %s: %d",szMes, nOffset_value);

		}
	}
	m_pUserMemory->m_iGoldenResult = 0;

	int iGoldenValueFlag = 0;
	int iGoldenSampleCount = 0;
	char szGoldenPath[260] = {0};
	char szTempFlag[32] = {0};

	wsprintf(szGoldenPath, "%s\\Config\\GoldenSample.ini", m_pPfMemory->GetRoot()); 
	m_pPfMemory->SetLog(nSite, "m_szGoldenPath: %s", szGoldenPath);
	iGoldenSampleCount = GetPrivateProfileInt("GoldenSampleSN", "Count", 0, szGoldenPath);
	iGoldenValueFlag = GetPrivateProfileInt("GoldenVerifyControl", "GoldenValue", 0, szGoldenPath);
	m_pPfMemory->SetLog(nSite, "GoldenValueFlag = %d", iGoldenValueFlag);
	if (1 == iGoldenValueFlag)
	{
		if (0 == (iGoldenSampleCount % 2))
		{
			AfxMessageBox("Golden Sample Count is Even !!! Can't Write T0 Value!!!");

			ExitProcess(0);
		}

		CleanGoldenT0AvgValue(nUsbID, nSite, iGoldenSampleCount, szGoldenPath);

		InitGoldenItem(nUsbID, nSite, iGoldenSampleCount, szGoldenPath);
	}


#ifndef _DEBUG
	char szCaseIniPath[MAX_PATH] = {0};
	wsprintf(szCaseIniPath, "%s\\Config\\Case.ini", m_pPfMemory->GetRoot()); 
	tstring strMD5 = _T("");
	tstring strMainPath = _T("");
	CString strTempPath = _T("");
	vector<tstring> vtFileType;
	vector<tstring> vtIgnore;

	strTempPath = m_pPfMemory->GetRoot();
	strMainPath = strTempPath;
	vtFileType.push_back(_T(".dll"));

	GetMD5(strMD5, strMainPath, vtFileType, vtIgnore);
	m_pPfMemory->SetLog(nSite, "strMD5: %s", strMD5.c_str());
	WritePrivateProfileString("SysConfig", "DllMD5", strMD5.c_str(), szCaseIniPath);

	OpenDolanDemo(nUsbID, nSite);

#endif

	GetSFCInfo(nUsbID, nSite);

	return 0;
}

BOOL CTest::InitGoldenItem(UINT nUsbID, UINT nSite, UINT iCount, const char* szConfigPath)
{
	int i = 0;
	int iIndex = 0;
	vector <ST_GoldenItem> vtGoldenItem;
	char *pBuffer = NULL;
	char *pPoint = NULL;
	ST_GoldenItem _GoldenItem;
	_GoldenItem.Init();
	char szTempGoldenItem[256] = {0};

	pBuffer = new char[65535];
	memset(pBuffer, NULL, sizeof(char)*65535);
	pPoint = pBuffer;

	GetPrivateProfileSection("Golden_Item", pBuffer, 65535, szConfigPath);
	pPoint = pBuffer;

	while(strlen(pPoint) != 0)
	{
		i = 0;
		while(pPoint[i] != '=')
		{
			i++;
		}
		pPoint[i] = '\0';

		_GoldenItem.strItemName = pPoint;
		vtGoldenItem.push_back(_GoldenItem);

		i = strlen(pPoint);
		i += 2;
		pPoint += i;
//		iIndex++;
	}

	for (i = 0; i < iCount; i++)
	{
		m_pUserMemory->m_vtGoldenItem.push_back(vtGoldenItem);

		for (int j = 0; j < vtGoldenItem.size(); j++)
		{
			m_pPfMemory->SetLog(nSite, "m_vtGoldenItem[%d].Name[%d]: %s", i, j , m_pUserMemory->m_vtGoldenItem[i][j].strItemName.c_str());
		}
	}

	delete [] pBuffer;
	return TRUE;
}

void CTest::WriteGoldenT0AvgValue(UINT nUsbID, UINT nSite, UINT iCount, const char* szConfigPath)
{
	int i = 0, j = 0;
	int iIndex = 0;
	char *pBuffer = NULL;
	char *pPoint = NULL;
	ST_GoldenItem _GoldenItem;
	_GoldenItem.Init();
	vector <ST_GoldenItem> vtGoldenT0Value;

	char szBarcode[18] = {0};
	char szTemp[32] = {0};
	char szTempValue[16] = {0};
	double dTempValue = 0.0;


	pBuffer = new char[65535];
	memset(pBuffer, NULL, sizeof(char)*65535);
	pPoint = pBuffer;

	GetPrivateProfileSection("Golden_Item", pBuffer, 65535, szConfigPath);
	pPoint = pBuffer;

	while(strlen(pPoint) != 0)
	{
		i = 0;
		while(pPoint[i] != '=')
		{
			i++;
		}
		pPoint[i] = '\0';

		_GoldenItem.strItemName = pPoint;
		vtGoldenT0Value.push_back(_GoldenItem);

		for (i = 0; i < iCount; ++i)
		{
			wsprintf(szTemp, "BarcodeGolden%d", i+1);

			GetPrivateProfileString("GoldenSampleSN", szTemp, NULL, szBarcode, sizeof(szBarcode), szConfigPath);

			wsprintf(szTemp, "%s_GoldenValue", szBarcode);

			GetPrivateProfileString(szTemp, pPoint, NULL, szTempValue, sizeof(szTempValue), szConfigPath);

			dTempValue = atof(szTempValue);

			vtGoldenT0Value[iIndex].vtGoldenData.push_back(dTempValue);

		}

		i = strlen(pPoint);
		i += 2;
		pPoint += i;
		iIndex++;
	}

	char szTempGoldenItem[32] = {0};
	char szTempGoldenValue[32] = {0};

	for_each(vtGoldenT0Value.begin(), vtGoldenT0Value.end(), calc_avg<ST_GoldenItem>());

	for (i = 0; i < iIndex; ++i)
	{
		sprintf(szTempGoldenItem, "%s", vtGoldenT0Value[i].strItemName.c_str());
		sprintf(szTempGoldenValue, "%.3f", vtGoldenT0Value[i].dAvgValue);

		WritePrivateProfileString("GoldenT0Value", szTempGoldenItem, szTempGoldenValue, szConfigPath);
	}

	return ;
}

void CTest::CleanGoldenT0AvgValue(UINT nUsbID, UINT nSite, const int iCount, const char* szConfigPath)
{
	int i = 0, j = 0;
	int iIndex = 0;
	char *pBuffer = NULL;
	char *pPoint = NULL;

	char szTempGoldenItem[32] = {0};
	char szTempGoldenValue[32] = {0};
	char szTempIndex[32] = {0};
	char szTempBarcode[32] = {0};
	char szTemp[64] = {0};

	pBuffer = new char[65535];
	memset(pBuffer, NULL, sizeof(char)*65535);
	pPoint = pBuffer;

	GetPrivateProfileSection("Golden_Item", pBuffer, 65535, szConfigPath);
	pPoint = pBuffer;

	while(strlen(pPoint) != 0)
	{
		i = 0;
		while(pPoint[i] != '=')
		{
			i++;
		}
		pPoint[i] = '\0';

		wsprintf(szTempGoldenItem, "%s", pPoint);
		wsprintf(szTempGoldenValue, " ");
		for (j = 1; j <= iCount; ++j)
		{
			sprintf(szTempIndex, "BarcodeGolden%d", j);
			GetPrivateProfileString("GoldenSampleSN", szTempIndex, NULL, szTempBarcode, sizeof(szTempBarcode), szConfigPath);

			sprintf(szTemp, "%s_GoldenValue",szTempBarcode);
			WritePrivateProfileString(szTemp, szTempGoldenItem, szTempGoldenValue, szConfigPath);
		}

		i = strlen(pPoint);
		i += 2;
		pPoint += i;
		iIndex++;
	}
}

void CTest::OpenDolanDemo(UINT nUsbID, UINT nSite)
{
	char szExePath[MAX_PATH] = {0};
	///220823 remove
	//sprintf(szExePath, "%s\\DolanDemo\\DolanDemo.exe", m_pPfMemory->GetRoot());
	//m_pPfMemory->SetLog(nSite, "ShellExecute: %s", szExePath);

	//ShellExecute(NULL, "open", szExePath, NULL, NULL, SW_HIDE);

	sprintf(szExePath, "%s\\MTCPTool.exe", m_pPfMemory->GetRoot());
	m_pPfMemory->SetLog(nSite, "ShellExecute: %s", szExePath);
	ShellExecute(NULL, "open", szExePath, NULL, NULL, SW_NORMAL);
}

BOOL CTest::CheckMachineisRetest(void)
{
	BOOL bRet = TRUE;
// 	uploaddll::datastr ^p = gcnew uploaddll::datastr();
// 	String ^strRet;
// 	char szRetFormMis[128] = {0};
// 
// 	strRet = p->checkmachineisretest();
// 
// 	strcpy(szRetFormMis, (char*)(void*)System::Runtime::InteropServices::Marshal::StringToHGlobalAnsi(strRet));
// 
// 	m_pPfMemory->SetLog(0, "checkmachineisretest RetFormMis: %s", szRetFormMis);
// 
// 	if (strcmp(szRetFormMis, "TRUE") != 0)
// 	{
// 		bRet = FALSE;
// 	}

	return bRet;
}

void CTest::GetSFCInfo(UINT nUsbID, UINT nSite)
{
	char szFilePath[MAX_PATH] = {0};
	sprintf(szFilePath, "%s\\Case\\RxSFCInfo.csv", m_pPfMemory->GetRoot());

	std::vector<std::string> vtSFCInfo;

	FILE *pFile = NULL;
	pFile = _tfopen(szFilePath, _T("rb"));

	if (NULL != pFile)
	{
		int iLine = 0;
		int iIndex = 0;
		char szLineData[1024] = {0};
		char *pszData = NULL;
		char *pszTemp = NULL;
		while (!feof(pFile))
		{
			if (NULL != fgets(szLineData, 1024, pFile))
			{
				pszData = szLineData;

				while (NULL != (pszTemp = strchr(pszData, ',')))
				{
					*pszTemp = 0;
					vtSFCInfo.push_back(pszData);
					pszData = pszTemp + 1;
				}
				vtSFCInfo.push_back(pszData);
			}
// 			for (int k = 1; k < vtSFCInfo.size(); ++k)
// 			{
// 				m_pUserMemory->m_mapSFCInfo[vtSFCInfo[0].c_str()].push_back(vtSFCInfo[k].c_str());
// 			}

			if (8 == vtSFCInfo.size())
			{
				for (int k = 1; k < 6; ++k)
				{
					m_pUserMemory->m_mapSFCInfo[vtSFCInfo[0].c_str()].push_back(vtSFCInfo[k].c_str());
				}

				for (int i = 0; i < vtSFCInfo[7].size(); ++i)
				{
					if (vtSFCInfo[7][i] == '\r' || vtSFCInfo[7][i] == '\n')
					{
						vtSFCInfo[7][i] = 0;
					}
				}

				m_pUserMemory->m_mapSFCInfo_LensSN[vtSFCInfo[6].c_str()] = vtSFCInfo[7].c_str();
			}

			vtSFCInfo.clear();
		}

		fclose(pFile);
	}
}