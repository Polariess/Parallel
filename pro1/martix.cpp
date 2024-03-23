#include <iostream>
#include <cstring>
#include <chrono>
#include <ratio>
using namespace std::chrono;
const int maxn=1e4+10; 
unsigned int sum[maxn];
unsigned int b[maxn][maxn];
unsigned int a[maxn];
unsigned int tb[maxn][maxn];
unsigned int ta[maxn];
void init(int n)
{
	memcpy(a,ta,n*4);
	for(int i=0;i<n;i++)
	{
		memcpy(b[i],tb[i],n*4);
	}
}
void mul_normal(int n)
{
	
	for(int i=0;i<n;i++)
	{
		sum[i]=0;
		for(int j=0;j<n;j++)
		{
			sum[i]+=b[j][i]*a[j];
		}
		
	}
//	unsigned int total=0;
//	for(int i=0;i<n;i++)
//	{
//		total+=sum[i];
//	}
//  return total; //��֤������ȷ��
}
void mul_pro(int n)
{
	for(int i=0;i<n;i++)
	{
		sum[i]=0;
	} 
	for(int j=0;j<n;j++)
	{
		for(int i=0;i<n;i++)
		{
			sum[i]+=b[j][i]*a[j];
		}
	}
//	unsigned int total=0;
//	for(int i=0;i<n;i++)
//	{
//		total+=sum[i];
//	}
//  return total; //��֤������ȷ��
}
void mul_zkpro_2(int n)
{
	//ѭ��չ��:�ܴ���Ϊ2^k,����Ҫ���б�Ե���� 
	for(int i=0;i<n;i++)
	{
		sum[i]=0;
	}
	for(int j=0;j<n;j++)
	{
		for(int i=0;i<n;i+=2)
		{
			sum[i]+=b[j][i]*a[j];
			sum[i+1]+=b[j][i+1]*a[j];
		} 
	}
//	unsigned int total=0;
//	for(int i=0;i<n;i++)
//	{
//		total+=sum[i];
//	}
//  return total; //��֤������ȷ��
}
void mul_zkpro_4(int n)
{
	//ѭ��չ��:�ܴ���Ϊ2^k,����Ҫ���б�Ե���� 
	for(int i=0;i<n;i++)
	{
		sum[i]=0;
	}
	for(int j=0;j<n;j++)
	{
		for(int i=0;i<n;i+=4)
		{
			sum[i]+=b[j][i]*a[j];
			sum[i+1]+=b[j][i+1]*a[j];
			sum[i+2]+=b[j][i+2]*a[j];
			sum[i+3]+=b[j][i+3]*a[j];
		} 
	}
//	unsigned int total=0;
//	for(int i=0;i<n;i++)
//	{
//		total+=sum[i];
//	}
//  return total; //��֤������ȷ��
}
void mul_zkpro_8(int n)
{
	//ѭ��չ��:�ܴ���Ϊ2^k,����Ҫ���б�Ե���� 
	for(int i=0;i<n;i++)
	{
		sum[i]=0;
	}
	for(int j=0;j<n;j++)
	{
		for(int i=0;i<n;i+=8)
		{
			sum[i]+=b[j][i]*a[j];
			sum[i+1]+=b[j][i+1]*a[j];
			sum[i+2]+=b[j][i+2]*a[j];
			sum[i+3]+=b[j][i+3]*a[j];
			sum[i+4]+=b[j][i+4]*a[j];
			sum[i+5]+=b[j][i+5]*a[j];
			sum[i+6]+=b[j][i+6]*a[j];
			sum[i+7]+=b[j][i+7]*a[j];
		} 
	}
//	unsigned int total=0;
//	for(int i=0;i<n;i++)
//	{
//		total+=sum[i];
//	}
//  return total; //��֤������ȷ��
}
int main()
{
	//n:���ݹ�ģ n*n:����������� �޶�ÿ�γ���ִ����������Ϊ10^9,TΪִ������ 
	int n=1<<12;
	int T=1000000000ll/(n*n);
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			tb[i][j]=i+j;
		}			
		ta[i]=i;
	}
	init(n);
	high_resolution_clock::time_point t1=high_resolution_clock::now();
	for(int i=1;i<=T;i++)
	{
	//ѡ��ͬ���Ż��㷨���ֱ��¼ƽ����ʱ 
	//	mul_normal(n);
		mul_pro(n);
	//	mul_zkpro_2(n);
	//	mul_zkpro_4(n);
	//	mul_zkpro_8(n); 
	}
	high_resolution_clock::time_point t2=high_resolution_clock::now();
	duration<double,std::milli> time_span=t2-t1;
	std::cout<<time_span.count()/T<<std::endl;
	return 0;
}
