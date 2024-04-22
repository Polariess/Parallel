#include<iostream>
#include<chrono>
#include<ratio>
#include<ctime>
#include<stdio.h>
#include<time.h>
#include <malloc.h>
#include<emmintrin.h>
#include<immintrin.h>
#include<cstring>
using namespace std;
const int N=100; 
const int MAX=1000;
float a[N][N];
void Print(float A[][N])
{
	for(int i=0;i<N;i++)
 	{
 		for(int j=0;j<N;j++)std::cout<<a[i][j]<<" ";
 		std::cout<<std::endl;
	}
	std::cout<<std::endl;	
}

void Prepare()
{
	srand (static_cast <unsigned> (time(0)));
	//����ԭʼ�����Ǿ���:Ҳ���㷨����������õ�����ȷ��� 
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<i;j++)
			a[i][j]=0.0;
		a[i][i]=1.0;
 		for(int j=i+1;j<N;j++)
 			a[i][j]=rand()%MAX;
 	}
	
 	//����ÿһ��ȫ������ǰ����У��ɿ����㷨���������
	//����������ľ���һ�����Ի�Ϊ��׼����������ʽ 
 	for(int k=0;k<N;k++)
 		for(int i=k+1;i<N;i++)
 			for(int j=0;j<N;j++)
 				a[i][j]=(int)(a[i][j]+a[k][j])%MAX;
}
void Solve_Normal(float A[][N])
{
	for(int k=0;k<N;k++)
		for(int i=k+1;i<N;i++)
		{
			float u=a[i][k]/a[k][k];
			for(int j=k+1;j<N;j++)
			//��ȷ����Ӧ����0 - n-1������֪��Ȼ��Ϊ�����Ǿ��󣬿��ں���ֱ�Ӵ��� 
				a[i][j]-=u*a[k][j];
			a[i][k] = 0.0; //��Ϊ�����Ǿ������½���0 
		}
}


void Solve_SSE_ALIGN(float A[][N]) {  //�����SSE�㷨
    for (int k = 0; k < N; k++) {
        __m128 t1 = _mm_set1_ps(A[k][k]);
        int j = k+1;

        //�ԷǶ��벿�ֽ��д��д��� 
        while ((long long)(&A[k][j])%16)
        {
            A[k][j] = A[k][j] / A[k][k];
            j++;
        }
        //cout << &m[k][j]<<endl;
        for ( ; j + 4 <= N; j += 4) {
            __m128 t2 = _mm_load_ps(&A[k][j]);   //�Ѷ��룬��load��storeָ��
            t2 = _mm_div_ps(t2, t1);
            _mm_store_ps(&A[k][j], t2);
        }
        //����ĩ���д��д��� 
        for (; j < N; j++) {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++) {
            __m128 vik = _mm_set1_ps(A[i][k]);
            j = k + 1;
            //����ͬ�� 
            while ((long long)(&A[k][j])%16)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
                j++;
            }
            for ( ; j + 4 <= N; j += 4) {
                __m128 vkj = _mm_load_ps(&A[k][j]);
                __m128 vij = _mm_loadu_ps(&A[i][j]);
                __m128 vx = _mm_mul_ps(vik, vkj);
                vij = _mm_sub_ps(vij, vx);
                _mm_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void Solve_SSE(float A[][N]) {
    for (int k = 0; k < N; k++) {
        __m128 t1 = _mm_set1_ps(A[k][k]);
        int j = 0;
        for (j = k + 1; j + 4 <= N; j += 4) {
            __m128 t2 = _mm_loadu_ps(&A[k][j]);
            t2 = _mm_div_ps(t2, t1);
            _mm_storeu_ps(&A[k][j], t2);
        }
        for (; j < N; j++) {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++) {
            __m128 vik = _mm_set1_ps(A[i][k]);
            for (j = k + 1; j + 4 <= N; j += 4) {
                __m128 vkj = _mm_loadu_ps(&A[k][j]);
                __m128 vij = _mm_loadu_ps(&A[i][j]);
                __m128 vx = _mm_mul_ps(vik, vkj);
                vij = _mm_sub_ps(vij, vx);
                _mm_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void Solve_AVX(float A[][N]) {
    for (int k = 0; k < N; k++) {
        __m256 t1 = _mm256_set1_ps(A[k][k]);
        int j = 0;
        for (j = k + 1; j + 8 <= N; j += 8) {
            __m256 t2 = _mm256_loadu_ps(&A[k][j]);
            t2 = _mm256_div_ps(t2, t1);
            _mm256_storeu_ps(&A[k][j], t2);
        }
        for (; j < N; j++) {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++) {
            __m256 vik = _mm256_set1_ps(A[i][k]);
            for (j = k + 1; j + 8 <= N; j += 8) {
                __m256 vkj = _mm256_loadu_ps(&A[k][j]);
                __m256 vij = _mm256_loadu_ps(&A[i][j]);
                __m256 vx = _mm256_mul_ps(vik, vkj);
                vij = _mm256_sub_ps(vij, vx);
                _mm256_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void Solve_AVX_ALIGN(float A[][N]) {
    for (int k = 0; k < N; k++) {
        __m256 t1 = _mm256_set1_ps(A[k][k]);
        int j = k+1;
        while ((long long)(&A[k][j]) % 32)
        {
            A[k][j] = A[k][j] / A[k][k];
            j++;
        }
        for ( ; j + 8 <= N; j += 8) {
            __m256 t2 = _mm256_load_ps(&A[k][j]);
            t2 = _mm256_div_ps(t2, t1);
            _mm256_store_ps(&A[k][j], t2);
        }
        for ( ; j < N; j++) {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++) {
            __m256 vik = _mm256_set1_ps(A[i][k]);
            j = k + 1;
            while ((long long)(&A[k][j]) % 32)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
                j++;
            }
            for ( ; j + 8 <= N; j += 8) {
                __m256 vkj = _mm256_load_ps(&A[k][j]);
                __m256 vij = _mm256_loadu_ps(&A[i][j]);
                __m256 vx = _mm256_mul_ps(vik, vkj);
                vij = _mm256_sub_ps(vij, vx);
                _mm256_storeu_ps(&A[i][j], vij);
            }
            for ( ; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void Solve_AVX_512(float A[][N]) {
    for (int k = 0; k < N; k++) {
        __m512 t1 = _mm512_set1_ps(A[k][k]);
        int j = 0;
        for (j = k + 1; j + 16 <= N; j += 16) {
            __m512 t2 = _mm512_loadu_ps(&A[k][j]);
            t2 = _mm512_div_ps(t2, t1);
            _mm512_storeu_ps(&A[k][j], t2);
        }
        for (; j < N; j++) {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++) {
            __m512 vik = _mm512_set1_ps(A[i][k]);
            for (j = k + 1; j + 16 <= N; j += 16) {
                __m512 vkj = _mm512_loadu_ps(&A[k][j]);
                __m512 vij = _mm512_loadu_ps(&A[i][j]);
                __m512 vx = _mm512_mul_ps(vik, vkj);
                vij = _mm512_sub_ps(vij, vx);
                _mm512_storeu_ps(&A[i][j], vij);
            }
            for (; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}

void Solve_AVX_512_ALIGN(float A[][N]) {
    for (int k = 0; k < N; k++) {
        __m512 t1 = _mm512_set1_ps(A[k][k]);
        int j = k+1;
        while ((long long)(&A[k][j]) % 64)
        {
            A[k][j] = A[k][j] / A[k][k];
            j++;
        }
        for ( ; j + 16 <= N; j += 16) {
            __m512 t2 = _mm512_load_ps(&A[k][j]);
            t2 = _mm512_div_ps(t2, t1);
            _mm512_store_ps(&A[k][j], t2);
        }
        for ( ; j < N; j++) {
            A[k][j] = A[k][j] / A[k][k];
        }
        A[k][k] = 1.0;
        for (int i = k + 1; i < N; i++) {
            __m512 vik = _mm512_set1_ps(A[i][k]);
            j = k + 1;
            while ((long long)(&A[k][j]) % 64)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
                j++;
            }
            for ( ; j + 16 <= N; j += 16) {
                __m512 vkj = _mm512_load_ps(&A[k][j]);
                __m512 vij = _mm512_loadu_ps(&A[i][j]);
                __m512 vx = _mm512_mul_ps(vik, vkj);
                vij = _mm512_sub_ps(vij, vx);
                _mm512_storeu_ps(&A[i][j], vij);
            }
            for ( ; j < N; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }
    }
}
int main()
{
	struct timespec start, end;
    long long elapsed_time; // �����洢������ʱ�䣨���룩

    // ��ӡ��ǰ���ݹ�ģ
    printf("��ǰ���ݹ�ģΪ:%d\n", N);

    // ƽ���㷨��ʱ
    Prepare();
    clock_gettime(CLOCK_MONOTONIC, &start);
    Solve_Normal(a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) * 1000LL + (end.tv_nsec - start.tv_nsec) / 1000000;
    printf("ƽ���㷨��ʱ:%lldms\n", elapsed_time);

    // SSE�޶����ʱ
    Prepare();
    clock_gettime(CLOCK_MONOTONIC, &start);
    Solve_SSE(a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) * 1000LL + (end.tv_nsec - start.tv_nsec) / 1000000;
    printf("SSE(�޶���)��ʱ:%lldms\n", elapsed_time);

    // SSE�����ʱ
    Prepare();
    clock_gettime(CLOCK_MONOTONIC, &start);
    Solve_SSE_ALIGN(a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) * 1000LL + (end.tv_nsec - start.tv_nsec) / 1000000;
    printf("SSE(����)��ʱ:%lldms\n", elapsed_time);
    
    // AVX�޶����ʱ
    Prepare();
    clock_gettime(CLOCK_MONOTONIC, &start);
    Solve_AVX(a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) * 1000LL + (end.tv_nsec - start.tv_nsec) / 1000000;
    printf("AVX(�޶���)��ʱ:%lldms\n", elapsed_time);

    // AVX�����ʱ
    Prepare();
    clock_gettime(CLOCK_MONOTONIC, &start);
    Solve_AVX_ALIGN(a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) * 1000LL + (end.tv_nsec - start.tv_nsec) / 1000000;
    printf("AVX(����)��ʱ:%lldms\n", elapsed_time);

//	return 0;
	
    // AVX-512�޶����ʱ
    Prepare();
    clock_gettime(CLOCK_MONOTONIC, &start);
    Solve_AVX_512(a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) * 1000LL + (end.tv_nsec - start.tv_nsec) / 1000000;
    printf("AVX-512(�޶���)��ʱ:%lldms\n", elapsed_time);

    // AVX-512�����ʱ
    Prepare();
    clock_gettime(CLOCK_MONOTONIC, &start);
    Solve_AVX_512_ALIGN(a);
    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) * 1000LL + (end.tv_nsec - start.tv_nsec) / 1000000;
    printf("AVX-512(����)��ʱ:%lldms\n", elapsed_time);
	return 0;
}
