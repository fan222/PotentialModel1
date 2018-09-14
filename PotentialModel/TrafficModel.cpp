#include "StdAfx.h"
#include "TrafficModel.h"



TrafficModel::TrafficModel(void)
{
}


TrafficModel::~TrafficModel(void)
{




}


  BOOL TrafficModel::TrafficSetting(SLane &LineInfo,SObstacle &ObsInfo,FrenetLaneMark &LaneMark)
  {
	  
	  if (LineInfo.IsOK==FALSE)//判断车道线是否有效
		  return FALSE;
	  LaneMark.LaneWidth=LineInfo.LaneWidth;

	  SLaneLine DKBaseLine;
	  memset(&DKBaseLine,0,sizeof(DKBaseLine));//frenet坐标系基准线

	  memset(&LaneMark,0,sizeof(FrenetLaneMark));
	

	  if (  LineInfo.LaneNo==3 || LineInfo.LaneNo==2)//挑选最右侧车道线为基准线
	  {
		  memcpy(&DKBaseLine,&LineInfo.LineR2,sizeof(LineInfo.LineR2));
		  if(LineInfo.LineR2.IsOK)
		  {
			  LaneMark.FrenetLineR2.LaneLineNo=1;
			  GetFrenetBaseLine(LineInfo.LineR2,LaneMark.FrenetLineR2);		      
			  FrenetLane[0]=LaneMark.FrenetLineR2;
		  } 
		  if(LineInfo.LineR1.IsOK)
		  {
			  GetFrenetLine(LineInfo.LineR1,DKBaseLine,LaneMark.FrenetLineR1,LaneMark.FrenetLineR2);
			  LaneMark.FrenetLineR1.LaneLineNo=2;
			  FrenetLane[1]=LaneMark.FrenetLineR1;
		  }
		  if (LineInfo.LineL1.IsOK)
		  {
			  GetFrenetLine(LineInfo.LineL1,DKBaseLine,LaneMark.FrenetLineL1,LaneMark.FrenetLineR2);
			  LaneMark.FrenetLineL1.LaneLineNo=3;
			  FrenetLane[2]=LaneMark.FrenetLineL1;
		  }
		  if (LineInfo.LineL2.IsOK)
		  {
			  GetFrenetLine(LineInfo.LineL2,DKBaseLine,LaneMark.FrenetLineL2,LaneMark.FrenetLineR2);
			  LaneMark.FrenetLineL2.LaneLineNo=4;
			  FrenetLane[3]=LaneMark.FrenetLineL2;
		  }
		 
	  }else if (LineInfo.LaneNo==1)
	  {
		  if (LineInfo.LineR1.IsOK)
		  {
			  memcpy(&DKBaseLine,&LineInfo.LineR1,sizeof(LineInfo.LineR1));
			  GetFrenetBaseLine(LineInfo.LineR1,LaneMark.FrenetLineR1);
			  LaneMark.FrenetLineR1.LaneLineNo=1;
			  FrenetLane[0]=LaneMark.FrenetLineR1;
		  }
		  if (LineInfo.LineL1.IsOK)
		  {
			  GetFrenetLine(LineInfo.LineL1,DKBaseLine,LaneMark.FrenetLineL1,LaneMark.FrenetLineR2);
			  LaneMark.FrenetLineL1.LaneLineNo=2;
			 FrenetLane[1]=LaneMark.FrenetLineL1;
		  }
		  if (LineInfo.LineL2.IsOK)
		  {
			  GetFrenetLine(LineInfo.LineL2,DKBaseLine,LaneMark.FrenetLineL2,LaneMark.FrenetLineR2);
			  LaneMark.FrenetLineL2.LaneLineNo=3;
			  FrenetLane[2]=LaneMark.FrenetLineL2;
		  }		
	  }

	  
		  GetFrenetObs(ObsInfo,DKBaseLine,LaneMark);

	 return TRUE;
  }

  BOOL TrafficModel::GetFrenetBaseLine(SLaneLine &line, FrenetLaneLine &FrenetLine)
  {
	  vector<fPoint> samples;//曲线拟合基本采样点
	  FrenetLine.IsOk = TRUE;

	  int LinePointNum = line.ValidNum;
	  FrenetLine.ValidNum=line.ValidNum;
	  double a[4]={0,0,0,0};
	  float s=0;
	  int cecl=(int)LinePointNum/5;
	  for (int i= 0; i<cecl*5; i=i+5)
	  {
		  FrenetLinePt tempPoint;

			for (int j=0;j<5;j++)
			{
				samples.push_back(line.LinePt[i+j]);

			}
			memset(&a, 0, sizeof(double)*4);
			Fitting(samples,a,4);
			samples.clear();

		
			for (int k=0;k<5;k++)
			{
				tempPoint.s = s;
				tempPoint.d = 0;		

				float dy = 3*a[3]*pow(line.LinePt[i+k].x,2) + 2*a[2]*line.LinePt[i+k].x + a[1];
				float d2y = 6*a[3]*line.LinePt[i+k].x + 2*a[2];
				if(i+k==0)
					tempPoint.kappa=0;
				else
					tempPoint.kappa = d2y/pow(1+dy*dy, (float)1.5);
				tempPoint.type = line.LaneLineClass;
				FrenetLine.FrenetLinePT.push_back(tempPoint);
				s += line.LinePt[i+k].Pt2PtDist(line.LinePt[i+k+1]);

				
				memset(&tempPoint,0,sizeof(tempPoint));
		   }
		 
	 }
		  
	  return TRUE;
  }
  BOOL TrafficModel::GetFrenetLine(SLaneLine &line,SLaneLine &DKBaseLine,FrenetLaneLine &FrenetLine,FrenetLaneLine &FrenetBaseLine)
  {
	  int Num = DKBaseLine.ValidNum;
	  int index = 0;
	  
	  float dtemp=1000;
	  FrenetLinePt temp;
	  FrenetLine.IsOk=TRUE;
	  FrenetLine.ValidNum=line.ValidNum;
	  for (int j=0;j<line.ValidNum;j++)
	  {
		  float d = 1000;
		  float s=0;
		  for(int i = 0; i<Num ; i++)
		  {
			  dtemp = line.LinePt[j].Pt2PtDist(DKBaseLine.LinePt[i]);
			  if (dtemp<=d)
			  {
				  d=dtemp;
				  index = i;
			  }else
				  break;
			  if(i!=0)
				  s+=DKBaseLine.LinePt[i-1].Pt2PtDist(DKBaseLine.LinePt[i]);
			  else 
				  continue;
			 
		 }
		  float  temps=0;
		  fPoint A=DKBaseLine.LinePt[index];
		  fPoint B=DKBaseLine.LinePt[index+1];
		  fPoint C=line.LinePt[j];
		  fPoint AB;
			     AB.x=B.x-A.x;
		         AB.y=B.y-A.y;
		  fPoint AC;
			     AC.x=C.x-A.x;
		         AC.y=C.y-A.y;
		  for (int i=0;i<10;i++)
		  {
			  fPoint insertpoint;
			         insertpoint.x=AB.x/10*i+A.x;
			         insertpoint.y=AB.y/10*i+A.y;
			  dtemp = insertpoint.Pt2PtDist(C);
			  if (dtemp<=d)
			  {
				  d=dtemp;
			  }
			  else
			  {
				  temps=sqrt(pow(AB.x*(i-1)/10,2)+pow(AB.y*(i-1)/10,2));
				  break;
			  }
		  }

		
		  temp.s=s+temps;
		  temp.d=d;
		  temp.kappa=(FrenetBaseLine.FrenetLinePT[index].kappa+FrenetBaseLine.FrenetLinePT[index].kappa)/2;
		  temp.type=line.LaneLineClass;
		  FrenetLine.FrenetLinePT.push_back(temp);
		  memset(&temp,0,sizeof(temp));
	  }
	  return TRUE;
  }

  BOOL TrafficModel::GetFrenetObs(SObstacle &Obs,SLaneLine &DKBaseLine,FrenetLaneMark &FrenetLine)
  {  

	  FrenetObsVeh      TempObsVeh;            //障碍车信息
	  FrenetObsHuman    TempObsHuman;         //行人
	  FrenetStaticObs    TempStaticObs;       //静态语义障碍物
	  int vehNum=0;
	  int humanNum=0;
	  int statcObsNum=0;

	  fPoint A;
	  fPoint B;
	  fPoint C;
	  fPoint AB;
	  fPoint AC;
	  fPoint insertpoint;

	  for (int j=0;j<Obs.ObstacleNum;j++)
	  {
		   float dtemp = 1000;
		  float d = 1000;
		  int BaseLinePtNum=DKBaseLine.ValidNum;
		  int index=0;
		  float s=0;
		  float theta=0;
		  for (int i=0;i<BaseLinePtNum;i++)
		  {
			  dtemp=Obs.Obs[j].Center.Pt2PtDist(DKBaseLine.LinePt[i]);
			  if (dtemp<=d)
			  {
				  d=dtemp;
				  index=i;
			  }
			  else
				  break;
			  if(i!=0)
				  s+=DKBaseLine.LinePt[i-1].Pt2PtDist(DKBaseLine.LinePt[i]);
			  else continue;
			
		  }


		  A=DKBaseLine.LinePt[index];

		  B=DKBaseLine.LinePt[index+1];
		
		  C=Obs.Obs[j].Center;
		  
		  AB.x=B.x-A.x;
		  AB.y=B.y-A.y;
		  
		  AC.x=C.x-A.x;
		  AC.y=C.y-A.y;
		  float temps=0;
		  for (int i=0;i<10;i++)
		  {
		
			  insertpoint.x=AB.x/10*i+A.x;
			  insertpoint.y=AB.y/10*i+A.y;
			  dtemp = insertpoint.Pt2PtDist(C);
			  if (dtemp<=d)
			  {
				  d=dtemp;
			  }
			  else
			  {
				  temps=sqrt(pow(AB.x*(i-1)/10,2)+pow(AB.y*(i-1)/10,2));
				  break;
			  }
		  }
		  s=s+temps;

		  if (Obs.Obs[j].ObstacleID==Vehicle)
		  {
			  TempObsVeh.s=s;
			  TempObsVeh.d=d;
			  TempObsVeh.v=Obs.Obs[j].v;
			  TempObsVeh.a=0;
			  TempObsVeh.len=5;
			  if (FrenetLine.FrenetLineR2.LaneLineNo==1)
			  {
				  TempObsVeh.theta=atan((1-FrenetLine.FrenetLineR2.FrenetLinePT[index].kappa*Obs.Obs[j].v)*tan(Obs.Obs[j].Heading));
			  }else
			  {
				  TempObsVeh.theta=atan((1-FrenetLine.FrenetLineR1.FrenetLinePT[index].kappa*Obs.Obs[j].v)*tan(Obs.Obs[j].Heading));
			  }
			  ObsVeh[vehNum]=TempObsVeh;
			  vehNum++;
		  }else if (Obs.Obs[j].ObstacleID == Human)
		  {
			  TempObsHuman.s=s;
			  TempObsHuman.d=d;
			  if (FrenetLine.FrenetLineR2.LaneLineNo==1)
			  {
				  TempObsVeh.theta=atan((1-FrenetLine.FrenetLineR2.FrenetLinePT[index].kappa*Obs.Obs[j].v)*tan(Obs.Obs[j].Heading));
			  }else
			  {
				  TempObsVeh.theta=atan((1-FrenetLine.FrenetLineR1.FrenetLinePT[index].kappa*Obs.Obs[j].v)*tan(Obs.Obs[j].Heading));
			  }
			  TempObsHuman.v=Obs.Obs[j].v;
			  TempObsHuman.a=0;
			  ObsHuman[humanNum]=TempObsHuman;
			  humanNum++;
		  }else
		  {
			  TempStaticObs.s=s;
			  TempStaticObs.d=d;
			  TempStaticObs.r=1.5;
			  TempStaticObs.sem_len=abs(Obs.Obs[j].LeftFront.Pt2PtDist(Obs.Obs[j].LeftBack));
			  StaticObs[statcObsNum]=TempStaticObs;
			  statcObsNum++;
		  }
	  }

	  
	  return TRUE;
  }


  BOOL TrafficModel::GetOptVelocity(SLane &LineInfo,SObstacle &ObsInfo, float t,float a)
  {




	  return TRUE;
  }
  BOOL TrafficModel::GetVelocityFile(float t,float a,vector<float> &velocity)
  {
	  return TRUE;
  }

  float TrafficModel::LaneMarkCost(float s,float d,FrenetLaneMark &LaneMark)
  {

	  return 0;
  }
  float TrafficModel::DynamicLaneRightCost(float s,float d,SObstacle &ObsInfo,FrenetLaneMark &LaneMark,AutoVeh &veh)
  {
	  return 0;
  }
  float TrafficModel::StaticLaneRightCost(float s,float d,SObstacle &ObsInfo,FrenetLaneMark &LaneMark,AutoVeh &veh)
  {
	  return 0;
  }
  float TrafficModel::HumanCost(float s,float d, SObstacle &ObsInfo)
  {
	  return 0;
  }

  float TrafficModel::GetStepDis(vector<float> &velocity,int LayerIndex)
  {




	  return 0;
  }
  BOOL TrafficModel::GetTreeMap(vector<vector< TreeNode >> &TreeMap)
  {
	  return TRUE;
  }
  BOOL TrafficModel::BellmanIteration(vector<vector< TreeNode >> &TreeMap)
  {
	  return TRUE;
  }
  BOOL TrafficModel::OptPolicy(vector<vector< TreeNode >> &TreeMap,vector<TreeNode> &OptNode)
  {
	  return TRUE;
  }

  // 曲线拟合
  void TrafficModel::Fitting(vector<fPoint> Samples, double* a, int nParam)
  {
	  int i, j, m, n;
	  m = Samples.size();
	  double Cond = 0;
	  n = nParam;
	  double A[N_PARAM][N_PARAM];
	  double b[N_PARAM];
	  for (i = 0; i < n; i ++) {
		  for (j = 0; j < n; j ++) {
			  A[i][j] = 0.;
		  }
		  b[i] = 0.;
	  }
	  A[0][0] = m;
	  for (i = 0; i < m; i ++) {
		  b[0] += Samples[i].y;
	  }
	  for (j = 1; j < n; j ++) {
		  for (i = 0; i < m; i ++) {
			  A[j][0] += pow(double(Samples[i].x), j);
			  A[n - 1][j] += pow(double(Samples[i].x), n + j - 1);
			  b[j] += pow(double(Samples[i].x), j) * Samples[i].y;
		  }
	  }
	  for (j = 1; j < n - 1; j ++) {
		  for (i = j; i < n - 1; i ++) {
			  A[i][j] = A[i + 1][j - 1];
		  }
	  }
	  for (i = 0; i < n - 1; i ++) {
		  for (j = i + 1; j < n; j ++) {
			  A[i][j] = A[j][i];
		  }
	  }
	  Matrix matrix(n);
	  matrix.m_SetMatrix(A);
	  Cond = matrix.m_TriangDecomp();
	  if (Cond < 0.) 
	  {
		  return;
	  }
	  matrix.m_SetVector(b);
	  matrix.m_BackSubstitute();
	  matrix.m_GetVector(b);
	  for (i = 0; i < n; i ++) {
		  a[i] = b[i];
	  }
	  return;
  }