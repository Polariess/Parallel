#include <iostream>
#include <cstring>
#include <chrono>
#include <ratio>
using namespace std::chrono;
const int maxn=1e7+10; 
unsigned int a[maxn];
unsigned int ta[maxn];
void init(int n)
{
	memcpy(a,ta,n*4);
}
void add_normal(int n)
{
	unsigned int sum=0;
	for(int i=0;i<n;i++)
	{
		sum+=a[i];	
	}
//  return sum; //��֤������ȷ��
}

void add_double_chain(int n)
{
	unsigned int sum=0,sum1=0,sum2=0;
	for(int i=0;i<n;i+=2)
	{
		sum1+=a[i];
		sum2+=a[i+1];
	}
	sum=sum1+sum2;
//  return sum; //��֤������ȷ��
}
void add_merge(int n)
{
	for(int m=n;m>1;m>>=1)
  	{
    	for(int i=0;i<m>>1;i++)
    	{
          	a[i]=a[i<<1]+a[(i<<1)+1];
    	}
  	}
  	//return a[0]; //��֤������ȷ��
}
int main()
{
	//n:���ݹ�ģ n:����������� �޶�ÿ�γ���ִ����������Ϊ10^9,TΪִ������
	int n=1<<17;
	//std::cout<<n<<std::endl;
	int T=1000000000ll/n;
	//std::cout<<T<<std::endl;
	for(int i=0;i<n;i++)
	{
		ta[i]=i+1;
	} 
	high_resolution_clock::time_point t1=high_resolution_clock::now();
	for(int i=1;i<=T;i++)
	{
		init(n);
		//ѡ��ͬ���Ż��㷨���ֱ��¼ƽ����ʱ 
		//add_normal(n);
		add_double_chain(n);
		//add_merge(n); 
	}
	high_resolution_clock::time_point t2=high_resolution_clock::now();
	duration<double,std::milli> time_span=t2-t1;
	std::cout<<time_span.count()/T<<std::endl;
	return 0;
}
