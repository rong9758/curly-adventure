#pragma once

#include "SpecDlg.h"
#include "Login.h"
#define NVMSIZE     1024
#define OTP_BufSize	49
#define MAX_PATH_LENGTH 520
/*
	Company			: HyVision System
	Product			: Platform HW	V5UFP MEMORY
	Description		: Inteface HW	V5UFP MEMORY
	Version			: 4.0.0.1
	
	Reporting date	: 2015.08.12
	Name			: Sung - Jin UK
	
	*
	* Do not modify.
	*
*/

class CTest
{

public:

	CHVS_InterFace_SW*					m_pSoftware;
	CHVS_InterFace_HW_V5FP*				m_pHardware;
	void*								m_pProtocol;
	CHVSPlatformMemory*					m_pPfMemory;

	CSpecDlg*							m_pSpecDlg;

	CUser_InterFace_Memory*				m_pUserMemory;
public:

	CTest();

	~CTest();

	BOOL OnInit(CHVS_InterFace_SW* hvs_SW, void* hvs_HW, void* hvs_Protocol, CHVSPlatformMemory* pMemory);

	UINT StartTest(UINT nUsbID, UINT nSite, UINT nTestNumber);
	CLogin* m_pLoginDlg;

	BOOL InitGoldenItem(UINT nUsbID, UINT nSite, UINT iCount, const char* szConfigPath);
	void WriteGoldenT0AvgValue(UINT nUsbID, UINT nSite, UINT iCount, const char* szConfigPath);
	void CleanGoldenT0AvgValue(UINT nUsbID, UINT nSite, const int iCount, const char* szConfigPath);

	void OpenDolanDemo(UINT nUsbID, UINT nSite);
	BOOL CheckMachineisRetest(void);

	void GetSFCInfo(UINT nUsbID, UINT nSite);
// 	char m_szConfigPath[MAX_PATH];
// 	CString m_szMainPath;
};

//extern //CTest m_Test;