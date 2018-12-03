
#include<fstream>
#include<iostream>
#include<iomanip>
#include<cstring>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<random> // has function for random no. generation.

using namespace std;

//--------------------------------------------
// SIMULATION PROGRAM
//-------------------------------------------
const int space_len=256; // length of lattice array
const int space_bre=256; // breadth of lattice array
const int NO_ENZYMES = 1;
const int TOT_TIMESTEPS = 1048576;
char space[space_len][space_bre];

void init_space()
{
    for(int i=0;i<space_len;i++)
        for(int j=0;j<space_bre;j++)
            space[i][j] = '#';
}

struct point
{
    int x; int y;
};


class rect
{
    int l,b;
    public:
        rect() {l=1;b=1;}
        rect(int len, int bre) {l=len; b=bre;}
        void putdim(int len, int bre) {l=len; b=bre;}
        void launch(point pos,char symbol);//pos represents coordinates of top left corner
};

void rect::launch(point pos, char symbol)// to put an obstacle on the lattice array
{
    for(int j=pos.y; j<(pos.y+b); j++)
        for(int i=0;i<1;i++)
            space[j][pos.x+i]=symbol;
}

int chk_empty(point position, int l, int b)//checks if the place where a new molecule has to be put is empty
{
    int indicator=1;
    for(int i=position.x;i<(position.x+l);i++)
        for(int j=position.y;j<(position.y+b);j++)
            if(space[j][i]!='#')
            {
                indicator = 0;
                break;
            }
	return indicator;
}

void release2_dex(float area_fraction, int NO_SUB, char symbol, int size_obj, point start_pos, int area_len, int area_bre, long *seed, float space_enz)// to release obstacles or enzymes into the lattice
{
    /* area fraction = [S] + [0] */
    if (area_fraction <= 0) return; // If no obstacles Simply returning;
    int l=size_obj;
    int b=size_obj;
    float no = (area_len*area_bre*area_fraction)/(l*b) -NO_SUB;// no. of rectangles(obstacles)which need to be put.
    if(no<=0) {cout<<"*Error*"<<endl<< "Too low area fraction"<<endl; return;}
    float p;// probability of choosing a location
    float freespace = (area_len*area_bre) - space_enz;
    p = no / freespace;
    cout<<"before";
    cout<<"**prob="<<p;
	cout<<"after";
    rect obj(l,b);
    int xoffset, yoffset;
    yoffset = (int)(rand()%(b-1));
    cout<<"\nyoffset"<<yoffset<<endl;

    int ctr=0;
    point pos;
    for(int yco=yoffset; yco+=(b); yco < (area_bre-b))
    {
        if(yco > (area_bre-b)) break;
        xoffset=(int)(rand()%(l-1));
        cout<<"\nxoffset"<<xoffset<<yco<<endl;
        for(int xco=xoffset; xco+=(l); xco < (area_len-1))
        {
            if(xco > (area_len-1)) break;
            cout<<xco;
            if((double)(rand()%1) < p) // i.e. probability is favourable
            {
                pos.x=xco;
                pos.y=yco;
                if(chk_empty(pos,l,b)==1)
                {
                    obj.launch(pos,symbol);
                    ctr++;
                }
            //if(ctr>=no) break;
            }
        }
    }
    cout<<"**ctr="<<ctr<<"**no="<<no<<endl;
}//end of release2_dex()

void point_enz_release(char symbol, int no, point start_pos, int area_len, int area_bre,long *seed)
{
    point launch_pos;
    for(int i=0;i<no;i++)
    {
        int xco,yco;
        do
        {
            xco=((int)(rand()%area_len))+start_pos.x; // chosing a set of coordinates
            yco=((int)(rand()%area_bre))+start_pos.y; // randomly
            launch_pos.x=xco;launch_pos.y=yco;
        }while(space[yco][xco]!='#');
        space[yco][xco]=symbol;
    }
}

class substrate
{
    public:
        point pos;//position of the particle on the lattice.
        point pos_real;//stores the actual position – coming due to cyclic boundary conditions
        static int pro_no;//stores the number of product molecules formed
        char state;
        point mother_enz;//store co-ordinate of that part of enz which made it product –appli. Only for products
        long time_birth; //stores the time when this product was formed – appli. only for products
        int res_time; // related to BINDING
        long n_jumps; // No. of actual jumps made for Products

        substrate(){state = 'S';mother_enz.x = -1;mother_enz.y = -1; time_birth = 0; res_time = 0; n_jumps =0;}
        void move(long *seed, float prob_rxn, int res_max, long cur_time);
        int chk_rxn(int xco, int yco, long *seed, float prob_rxn, long cur_time)
        {
            if(state =='P') return 0;
            if((space[yco][xco]=='E')&& ((double)(rand()%1)<=prob_rxn))
            {
                pro_no++;state='P';
                mother_enz.x=xco; mother_enz.y=yco; time_birth = cur_time;
                n_jumps = 0;
                return 1;
            }
            else return 0;
        }
};

int substrate::pro_no;

void substrate::move(long *seed,float prob_rxn,int res_max,long cur_time)
{
	if (res_time>0)
	{
		res_time--;
		return;
	}
	int m;
	//generate a random number and store it in m(an integer variable)
	m = (int)(rand()%4+1);

	switch(m)
	{
		//UP//
		case 1:if (pos.y==0)	//taking care of of boundary
		{
			if((state!='P')&&(chk_rxn(pos.x,space_bre-1,seed,prob_rxn,cur_time)==1))
			{
				space[pos.y][pos.x]=state;
			}
			else if ((space[space_bre-1][pos.x])=='#')
			{
				space[pos.y][pos.x]='#';
				pos.y=space_bre-1;
				space[pos.y][pos.x]=state;
				(pos_real.y)--;
				if (state=='P')
				{
					n_jumps++;
				}
			}
			else if((space[space_bre-1][pos.x])=='O')
			{
				res_time = res_max;
			}
		}
		else
		{
			if ((state!='P')&&(chk_rxn(pos.x,pos.y-1,seed,prob_rxn,cur_time)==1))
			{
				space[pos.y][pos.x]=state;
			}
			else if ((space[pos.y-1][pos.x])=='#')
			{
				space[pos.y][pos.x]='#';
				(pos.y)--;
				space[pos.y][pos.x]=state;
				(pos_real.y)--;
				if (state=='P')
				{
					n_jumps++;
				}
			}

			else if ((space[pos.y-1][pos.x])=='O')
			{
				res_time=res_max;
			}
		}
		break;

		//DOWN//
		case 2:if (pos.y==(space_bre-1))
		{
			if ((state!='P')&&(chk_rxn(pos.x,0,seed,prob_rxn,cur_time)==1))
			{
				space[pos.y][pos.x]=state;
			}
			else if ((space[0][pos.x])=='#')
			{
				space[pos.y][pos.x]='#';
				pos.y=0;
				space[pos.y][pos.x]=state;
				(pos_real.y)++;
				if (state=='P')
				{
					n_jumps++;
				}
			}
			else if((space[0][pos.x])=='O')
			{
				res_time = res_max;
			}
		}
		else
		{
			if ((state!='P')&&(chk_rxn(pos.x,pos.y+1,seed,prob_rxn,cur_time)==1))
			{
				space[pos.y][pos.x]=state;
			}
			else if ((space[pos.y+1][pos.x])=='#')
			{
				space[pos.y][pos.x]='#';
				(pos.y)++;
				space[pos.y][pos.x]=state;
				(pos_real.y)++;
				if (state=='P')
				{
					n_jumps++;
				}
			}

			else if ((space[pos.y+1][pos.x])=='O')
			{
				res_time=res_max;
			}
		}
		break;

		//RIGHT//
		case 3:if (pos.x==(space_len-1))
		{
			if ((state!='P')&&(chk_rxn(0,pos.y,seed,prob_rxn,cur_time)==1))
			{
				space[pos.y][pos.x]=state;
			}
			else if ((space[pos.y][0])=='#')
			{
				space[pos.y][pos.x]='#';
				pos.x=0;
				space[pos.y][pos.x]=state;
				(pos_real.x)++;
				if (state=='P')
				{
					n_jumps++;
				}
			}
			else if((space[pos.y][0])=='O')
			{
				res_time = res_max;
			}
		}
		else
		{
			if ((state!='P')&&(chk_rxn(pos.x+1,pos.y,seed,prob_rxn,cur_time)==1))
			{
				space[pos.y][pos.x]=state;
			}
			else if ((space[pos.y][pos.x+1])=='#')
			{
				space[pos.y][pos.x]='#';
				(pos.x)++;
				space[pos.y][pos.x]=state;
				(pos_real.x)++;
				if (state=='P')
				{
					n_jumps++;
				}
			}

			else if ((space[pos.y][pos.x+1])=='O')
			{
				res_time=res_max;
			}
		}
		break;

		//LEFT//
		case 4:if (pos.x==0)
		{
			if ((state!='P')&&(chk_rxn(space_len-1,pos.y,seed,prob_rxn,cur_time)==1))
			{
				space[pos.y][pos.x]=state;
			}
			else if ((space[pos.y][space_len-1])=='#')
			{
				space[pos.y][pos.x]='#';
				pos.x=space_len-1;
				space[pos.y][pos.x]=state;
				(pos_real.x)--;
				if (state=='P')
				{
					n_jumps++;
				}
			}
			else if((space[pos.y][space_len-1])=='O')
			{
				res_time = res_max;
			}
		}
		else
		{
			if ((state!='P')&&(chk_rxn(pos.x-1,pos.y,seed,prob_rxn,cur_time)==1))
			{
				space[pos.y][pos.x]=state;
			}
			else if ((space[pos.y][pos.x-1])=='#')
			{
				space[pos.y][pos.x]='#';
				(pos.x)--;
				space[pos.y][pos.x]=state;
				(pos_real.x)--;
				if (state=='P')
				{
					n_jumps++;
				}
			}

			else if ((space[pos.y][pos.x-1])=='O')
			{
				res_time=res_max;
			}
		}
		break;

	}	//end of switch
}		//end of move()

void sub_release(substrate sub[],char symbol,int start_index,int no,point start_pos,int area_len,int area_bre,long *seed)
{
	point launch_pos;
	for (int i = start_index; i < no; i++)
	{
		int xco,yco;
		do
		{
			xco = ((int)(rand()%area_len))+start_pos.x;		//chosing a set of co-ordinates randomly
			yco = ((int)(rand()%area_len))+start_pos.y;
			launch_pos.x = xco;
			launch_pos.y = yco;
		}while(space[yco][xco]!='#');
		sub[i].pos = launch_pos;
		sub[i].pos_real = launch_pos;
		space[yco][xco] = symbol;
	}
}

double cal_avg_jumps(substrate sub[], int no, long cur_time)
{
	double sum =0;
	int no_pro =0;
	for(int i =0; i < no; i++)
	{
		if(sub[i].state == 'P')
		{
			//cout<<cur_time<<” “<<sub[i].time_birth<<” “<<sub[i].n_jumps<<’#’;
			float age = cur_time - sub[i].time_birth;
			//cout<<sub[i].n_jumps/age<<”*”;
			if(age>0)
			{
				sum += (sub[i].n_jumps/age);
				no_pro++;
			}
		}
	}//end of for
	if(no_pro>0)
	{
		double avg_jumps = sum/no_pro;
		//cout<<endl;
		return avg_jumps;
	}
	else return 0;
}

double cal_diff_length(substrate sub[], int no)
{
	double sum =0;
	int no_pro =0;
	for(int i =0; i < no; i++)
	{
		if(sub[i].state == 'P')
		{
			double r_square =pow((sub[i].pos_real.x - sub[i].mother_enz.x),2.0) + pow((sub[i].pos_real.y  - sub[i].mother_enz.y),2.0);
			sum += r_square;
			no_pro++;
		}
	}//end of for
	if(no_pro>0)
	{
		double diff_length = sqrt(sum/(double)no_pro);
		return diff_length;
	}
	else return 0;
}

double cal_rad_gyration(substrate sub[], int no)
{
	double sum =0;
	int no_pro =0;
	//calculating center of mass
	float xcm = 0, ycm = 0;
	for(int i =0; i < no; i++)
	{
		if(sub[i].state == 'P')
		{
			xcm += sub[i].pos_real.x;
			ycm += sub[i].pos_real.y;
			no_pro++;
		}
	}//end of for
	if(no_pro == 0) return 0;
	xcm = xcm/no_pro; ycm=ycm/no_pro;
	//cout<<"## "<<xcm<<" "<<ycm<<endl;
	for(int i = 0; i<no;i++)
	{
	if(sub[i].state == 'P')
		{
			double r_square = (sub[i].pos_real.x - xcm)*(sub[i].pos_real.x - xcm) + (sub[i].pos_real.y - ycm)*(sub[i].pos_real.y - ycm);
			sum += r_square;
		}
	}//end of for
	if(no_pro>0)
	{
		double rad_gyration = sqrt(sum/(double)no_pro);
		return rad_gyration;
	}
	else return 0;
}

double QSP(long results[])
{
	int NS = 0, NP = 0;
	int NSP = 0,NSS = 0, NPP = 0;
	for(int i=0; i<(space_len -1);i++)
	{
		for(int j=0; j<(space_bre -1);j++)
		{
			if(space[i][j] == 'S')
			{
				NS++;
				if(space[i][j] == space[i][j+1]) NSS++;
				if(space[i][j] == space[i+1][j]) NSS++;
				if(space[i][j+1] == 'P') NSP++;
 				if(space[i+1][j] == 'P') NSP++;
			}
			if(space[i][j] == 'P')
			{
				NP++;
				if(space[i][j] == space[i][j+1]) NPP++;
				if(space[i][j] == space[i+1][j]) NPP++;
				if(space[i][j+1] == 'S') NSP++;
 				if(space[i+1][j] == 'P') NSP++;
			}
		}
		if(space[i][space_bre -1]=='S') NS++;
		if(space[i][space_bre -1]=='P') NP++;
	}

	for(int j=0; j<space_len;j++)
	{
		if(space[space_len -1][j]=='S') NS++;
		if(space[space_len -1][j]=='P') NP++;
	}
	results[0] = NSS;
	results[1] = NSP;
	results[2] = NPP;
	double Qsp_val = -1;
	if(NSP>0)
	{
		float fac = (float)2*NS*NP/(float)(NS*NS +NP*NP);
		Qsp_val = ((NSS+NPP)/(float)NSP)*fac;
	}
	 return Qsp_val;
}

void filewrite(char fname[], int no_entries)//for writing reaction files on to the disk
{
	ofstream fout;
	fout.open(fname);
	int i;
	fout<<no_entries<<'\n';//output
	for(i=0;i<space_len;i++)
	{
		for(int j=0; j<space_bre;j++)
		{
			if((space[i][j]!='#')&&(space[i][j]!='O')) fout<<j<<" "<<i<<" "<<space[i][j]<<"\n";
		}

	}
	fout.close();
}

void run_sim(substrate sub[], int start_index, int no, long timesteps, long *seed, float prob_rxn, int no_entries, int res_time, char fn_1g_out[], char fn_out[])
{
	char fname[][20] = {"posout1.txt", "posout2.txt", "posout3.txt", "posout4.txt", "posout5.txt", "posout6.txt", "posout7.txt", "posout8.txt", "posout9.txt", "posout10.txt", "posout11.txt", "posout12.txt", "posout13.txt", "posout14.txt", "posout15.txt", "posout16.txt", "posout17.txt", "posout18.txt", "posout19.txt", "posout20.txt"};
	int k=0;
	int log_out_indicator = 1;
	int out_indicator =10, out_ctr =0;
	ofstream fout_1g, fout_out;
	fout_1g.open(fn_1g_out);
	fout_out.open(fn_out);
	long results[3];
	for(int j=1;j<=timesteps;j++)
	{
		for(int i = start_index;i<(start_index+no);i++)
		{
			sub[i].move(seed,prob_rxn,res_time,j);
			if((j==log_out_indicator)&&(i==(start_index +no -1)))
			{
				fout_1g<<j<<"\t"<<sub[i].pro_no<<"\t"<<cal_rad_gyration(sub,no)<<"\t"<<cal_avg_jumps(sub,no,j)<<"\t";
				float Qsp_val = QSP(results);
				fout_1g<<results[0]<<"\t"<<results[1]<<"\t"<<results[2]<<"\t"<<Qsp_val<<endl;
				log_out_indicator*=2;
			}
			if(((j%log_out_indicator)==0)&&(i==(start_index +no -1)))
			{
				fout_out<<j<<"\t"<<sub[i].pro_no<<"\t"<<cal_rad_gyration(sub,no)<<"\t"<<cal_avg_jumps(sub,no,j)<<"\t";
				float Qsp_val = QSP(results);
				fout_out<<results[0]<<"\t"<<results[1]<<"\t"<<results[2]<<"\t"<<Qsp_val<<endl;
				out_ctr++;
				if((out_ctr%10)==9)
				{
					out_indicator*= 10; out_ctr++;
				}
			}
		}
	}
	fout_1g.close();
	fout_out.close();
}//end of run_sim

main()
{
	//inputting the required data
	float area_fraction ,sub_area_fraction,no_dex,prob_rxn;
	int size_dex;
	int res_time = 0;
	const int virtual_no_enz = NO_ENZYMES;
	int no_ens;
	cout<<"Enter the value of residence time"<<endl;
	cin >>res_time;
	cout<<"Enter the value of area_fraction"<<endl;
	cin>>area_fraction;
	cout<<"Enter the value of sub_area_fraction"<<endl;
	cin>>sub_area_fraction;
	cout<<"Enter the value of number of enzymes"<<endl;
	cin>>no_ens;
	//inputiing seed for random number generation
	long val; 	//seed
	cout<<"Enter the seed for random number"<<endl;
	cin>>val;	//Entering the random number seed

	//inputting part ends here
	size_dex=1,prob_rxn=1;
	const int virtual_no_sub = (int)(sub_area_fraction*space_len*space_bre);

	//Formatted output

	cout<<"Lattice_Length = "<<space_len<<"Lattice_Breadth ="<<space_bre<<"OBS :"<<size_dex<<"units "<<endl;
	cout<<"Resident Time ="<<res_time<<endl;
	cout<<"SUBSTRATE: "<<"1 units "<<virtual_no_sub<<"Molecules"<<endl;
	cout<<"ENZYME: "<<"1 units "<<(int)virtual_no_enz<<"Molecules"<<endl;
	cout<<"PROBABILITY: "<<prob_rxn<<endl;
	cout<<"RANDOM NO. SEED:"<<val<<endl;
	cout<<"AREA FRACTION OCCUPIED:"<<area_fraction<<endl;

	for (int ens = 1; ens <= no_ens; ens++)
	{
		cout<<"Ensemble#"<<ens<<"seed="<<val<<endl;
		substrate sub[10000];
		substrate :: pro_no = 0;
		init_space();	//initializing the lattice with '#'
		point po;
		po.x=0;
		po.y=0;
		space[(space_len)/2 - 1][(space_bre)/2 - 1] = 'E'; 		//Putting enzyme at the centre of the lattice
		//point_enz_release('E',virtual_no_enz,po,space_len,space_bre,&val);

		release2_dex(area_fraction,virtual_no_sub,'O',size_dex,po,space_len,space_bre,&val,(float)virtual_no_enz);

		//releasing substrate molecules
		sub_release(sub,'S',0,virtual_no_sub,po,space_len,space_bre,&val);
		int no_entries =0;
		//counting no of entries for reaction files
		for (int i = 0; i < space_len; i++)
		{
			for (int j = 0; i < space_bre; j++)
			{
				if ((space[i][j]!='#')&&(space[i][j]!='O'))
				{
					no_entries++;
				}
			}
		}

		char name_lg[20]= "log_out";
		char name_out[20] = "out";
		char num[10];
		sprintf(num,"%d",ens);
		strcat(name_lg,num);
		strcat(name_lg,".txt");
		strcat(name_out,num);
		strcat(name_out,".txt");

		//****************RUNNING THE SIMULATION **********************//
		run_sim(sub,0,virtual_no_sub,TOT_TIMESTEPS,&val,prob_rxn,no_entries,res_time,name_lg,name_out);

	}
}
