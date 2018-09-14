#pragma once
#include "stdafx.h"
#include "struct_fun.h"
#include "Matrix.h"
#define LayerNum 10;
#define BranchNum 3;
#define dt 0.1;





struct TreeNode
{
	float s;
	float d;
	float theta;
	float v;
	float cost;
	TreeNode()
	{
		s=0;
		d=0;
		theta=0;
		v=0;
		cost=0;
	}
};

struct FrenetLinePt 
{
	float s;
	float d;
	float kappa;
	int type;
	FrenetLinePt()
	{
		s=0;
		d=0;
		kappa=0;
		type=-1;
	}
};
struct FrenetLaneLine
{
	BOOL IsOk;
	int LaneLineNo;
	int ValidNum;
	vector<FrenetLinePt> FrenetLinePT;
	FrenetLaneLine()
	{
		IsOk =FALSE;
		FrenetLinePT.resize(600);
	}
};
struct FrenetObsVeh
{
	float s;
	float d;
	float theta;
	float v;
	float a;
	float len;
	FrenetObsVeh()
	{
		s=0;
		d=0;
		theta=0;
		v=0;
		a=0;
		len=0;
	}
};
struct FrenetObsHuman
{
	float s;
	float d;
	float theta;
	float v;
	float a;
	FrenetObsHuman()
	{
		s=0;
		d=0;
		theta=0;
		v=0;
		a=0;
	}
};
struct FrenetStaticObs
{
	float s;
	float d;
	float r;
	float sem_len;
	FrenetStaticObs()
	{
		s=0;
		d=0;
		r=0;
		sem_len=0;
	}

};

struct FrenetLaneMark
{
	BOOL IsOk;
	int LaneNo;//��ǰ�������ھ��Գ����ţ������ҳ���Ϊ0�����������ۼ�
	double LaneWidth;//������
	FrenetLaneLine		FrenetLineR2;			//����ο����Ҳ�ڶ���
	FrenetLaneLine		FrenetLineR1;			//����ο����Ҳ��һ��
	FrenetLaneLine		FrenetLineL1;			//����ο�������һ��
	FrenetLaneLine		FrenetLineL2;			//����ο������ڶ���
	


};

class TrafficModel
{
public:
	AutoVeh veh;
	//vector<vector< TreeNode >> TreeMap;//������������������Map,Layer��Node��
	//vector<float> velocity;//
	FrenetLaneLine   FrenetLane[4];
	FrenetObsVeh     ObsVeh[20];            //�ϰ�����Ϣ
	FrenetObsHuman   ObsHuman[20];         //����
	FrenetStaticObs  StaticObs[20];       //��̬�����ϰ���

	TrafficModel(void);
	~TrafficModel(void);


	BOOL GetOptVelocity(SLane &LineInfo,SObstacle &ObsInfo, float t,float a);//��ȡ����Ŀ����ٶȼ�ʱ��
	BOOL GetVelocityFile(float t,float a,vector<float> &velocity);//ʹ�ó���״̬����ģ��Ԥ�⳵���ٶ�


    BOOL TrafficSetting(SLane &LineInfo,SObstacle &ObsInfo,FrenetLaneMark &tLaneMark);//����������Ϣ״̬�趨�����ѿ���ת��FreNet��
	BOOL GetFrenetBaseLine(SLaneLine &line, FrenetLaneLine &FrenetLine);
	BOOL GetFrenetLine(SLaneLine &line,SLaneLine &DKBaseLine,FrenetLaneLine &FrenetLine,FrenetLaneLine &FreneBasetLine);
	BOOL GetFrenetObs(SObstacle &Obs,SLaneLine &DKBaseLine,FrenetLaneMark &FrenetLine);





	float LaneMarkCost(float s,float d,FrenetLaneMark &LaneMark);//�������Ƴ�����
	float DynamicLaneRightCost(float s,float d,SObstacle &ObsInfo,FrenetLaneMark &LaneMark,AutoVeh &veh);//��̬����ռ���Ƴ�
	float StaticLaneRightCost(float s,float d,SObstacle &ObsInfo,FrenetLaneMark &LaneMark,AutoVeh &veh);//��̬����ռ���Ƴ�
	float HumanCost(float s,float d, SObstacle &ObsInfo);//����ռ���Ƴ�

	float GetStepDis(vector<float> &velocity,int LayerIndex);//��ȡĳһ��Ĳ���
	BOOL GetTreeMap(vector<vector< TreeNode >> &TreeMap);//����״̬��
	BOOL BellmanIteration(vector<vector< TreeNode >> &TreeMap);//������ֵ����
	BOOL OptPolicy(vector<vector< TreeNode >> &TreeMap,vector<TreeNode> &OptNode);//����·��������

	void Fitting(vector<fPoint> Samples, double* a, int nParam);//�������
};


