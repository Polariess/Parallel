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
//  return sum; //验证程序正确性
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
//  return sum; //验证程序正确性
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
  	//return a[0]; //验证程序正确性
}
int main()
{
	//n:数据规模 n:单轮运算次数 限定每次程序执行总运算量为10^9,T为执行轮数
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
		//选择不同的优化算法，分别记录平均用时 
		//add_normal(n);
		add_double_chain(n);
		//add_merge(n); 
	}
	high_resolution_clock::time_point t2=high_resolution_clock::now();
	duration<double,std::milli> time_span=t2-t1;
	std::cout<<time_span.count()/T<<std::endl;
	return 0;
}
