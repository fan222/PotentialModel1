#pragma once
#include "stdafx.h"


#define	MaxSpeed		18

#define LinePtNumMax	600		// �������ߵ�������
#define CARLENGTH 3.8	//����
#define CARHALFWIDTH 1.1	//����
#define Step	0.25   // �ֱ��ʣ���λm
#define ObstacleMaxNum	100		// �ϰ���������
#define N_PARAM		4		//	�ֶ������߶���ʽ����		



// ��������ϵ�еĵ�����
class fPoint
{
public:
	float x;
	float y;
	fPoint(){
	}
	~fPoint(){
	}
	fPoint(float a,float b){
		x = a;y = b ;
	}
	float Pt2PtDist(fPoint &tmp)
	{
		return sqrt(pow((x-tmp.x),2)+pow((y-tmp.y),2));
	}
};

// ������
class SLaneLine
{
public:
	BOOL			IsOK;
	int				LaneLineNo;
	int				LaneLineClass;			// 0 - �ɴ�Խ�� 1 - ���ɴ�Խ�� 3 - ������
	int				ValidNum;
	fPoint			LinePt[LinePtNumMax];
	//vector<fPoint>      LinePt;
	~SLaneLine(){}
	SLaneLine()
	{
		IsOK = FALSE;
		LaneLineNo = 9;
		LaneLineClass = 0;
		ValidNum = 0;
	}


};


// �ϰ�������
enum EObsMotionType
{	
	Static = 0,	   // static. 
	Moving = 1,	   // moving.
	TooLow = -5,     // the height difference is too small. 
	Unkown = -4,     // unkown motion type. 
	FewPts = -3,     // the obstacle contains few points. 
	UnMatch = -2,    // the current obstacle cannot find corresponding points in previous frame within some meters. 
	BadICP = -1      // icp results is low quality.
};

enum EObstacleClass
{
	// 0 un_know 1 vehicle 2 Motorcycle 3 Truck 4 Pedestrian 9 Bicycle.
	unknow = 0,
	Vehicle = 1, 
	Motorcycle = 2, 
	Truck = 3,       
	Human = 4,  
	warning = 5,	//	��ʾ��
	Cone = 6,	//	׶�ͱ�	   
	Bicycle = 9,  	
};

// �����ṹ
class SLane
{
public:
	BOOL			IsOK;			//TRUE-���ã�FALSE-������
	BOOL            IsVirtual;		//
	BOOL            IsSplit;
	int             CurRoadPriority;//
	int             LaneNo;	    	//��ǰ�������ھ��Գ����ţ������ҳ���Ϊ0�����������ۼ�
	double			LaneWidth;		//�������
	double          SplitDis;
	SLaneLine		LineR2;			//����ο����Ҳ�ڶ���
	SLaneLine		LineR1;			//����ο����Ҳ��һ��
	SLaneLine		LineL1;			//����ο�������һ��
	SLaneLine		LineL2;			//����ο������ڶ���
	SLaneLine       GlobalGuideLine;//ȫ��������
	double			LaneMaxSpeed;	//��ǰ֡������ʻ������ٶ�
	SLane(){
		IsSplit = FALSE;
		IsOK = FALSE;
		IsVirtual = FALSE;
		CurRoadPriority = 0;
		LaneNo = 9;
		LaneWidth = 0;
		SplitDis = 1000.0;
		LaneMaxSpeed = MaxSpeed;
	}

};

struct SObstacleInfo
{
	EObsMotionType	 MotionState;	//	�ϰ���״̬()
	fPoint 			 Center;		//	�ϰ���BOX�����ĵ��ھֲ�����ϵ�е�����
	fPoint 		     LeftBack;		//	�ϰ���BOX�������ھֲ�����ϵ�е�����
	fPoint 		     LeftFront;		//	�ϰ���BOX����ǰ���ھֲ�����ϵ�е�����
	fPoint 		     RightBack;		//	�ϰ���BOX���Һ���ھֲ�����ϵ�е�����
	fPoint 		     RightFront;	//	�ϰ���BOX����ǰ���ھֲ�����ϵ�е�����
	float			 Heading;		//	�ֲ�����ϵ�е��ٶȳ���,��λ�Ƕ�
	float			 v;				//	�ֲ�����ϵ�е�����ƶ��ٶ�
	float			 a;				//	�ֲ�����ϵ�еľ����ƶ����ٶȣ��޷��ṩ��0 
	int				 ObsLaneNo;		//  �ϰ������ڳ�����Ϣ
	EObstacleClass	 ObstacleID;	//	�ϰ��������Ϣ
	union UAddInfo						//  �������Ϣ����һλΪnTrackID;   ����λ�ñ����� 
	{
		struct 
		{
			int nTrackID;
			int nAge;
			int nTest; 
			int nObsTypeAge;  // �ϰ��ﱻ�۲⵽�������ԵĴ�����Ŀǰ��ֵ����IFV��LUX�ں�ʱʹ��
			int nObsSource;	  // 1:FLUX 2:RLUX 3:HDL 4:IFV 5:LUTM 6:RUTM 7:OutofRoadEdge
			int nReservedArr[5];        // reserved 
		};
		int data[10];
		union UAddInfo()
		{
			nAge = 0;
			nTrackID = -1;
			nObsTypeAge = 0;
		}
	}AddInfo;

};

// �ϰ���ṹ
// ObstacleMaxNum����ϰ�������
struct SObstacle
{
	BOOL			IsOK;						//  TRUE �C �� FALSE �C û��
	int 			ObstacleNum;				//	�ϰ�������
	SObstacleInfo	Obs[ObstacleMaxNum];		//	�ϰ�������
};

//����״̬�ṹ��
struct AutoVeh
{
	float s;
	float d;
	float theta;
	float v;
};
