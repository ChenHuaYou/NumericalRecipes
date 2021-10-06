
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(float **a, int n, float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	float big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

#define NR__gaussj(T)\
void NR__##T##__gaussj(NR__Matrix__##T *a, NR__Matrix__##T *b){\
    NR__Matrix__##T *I = NR__Matrix__##T##__eye(a->nrow);\
    NR__Matrix__##T *A = NR__Matrix__##T##__cat_multi_matrix(3,a,b,I);\
    for(int i=0; i<A->ncol; i++){\
        int idx = NR__Matrix__##T##__max_along_col(A, i);\
        NR__Matrix__##T##__row_swap(A, idx, i);\
        NR__Matrix__##T##__row_linear_combination(A, i, 1/A->data[i][i], i, 0);\
        for(int j=0; j<A->nrow; j++){\
            if(j==i) continue;\
            NR__Matrix__##T##__row_linear_combination(A, j, 1, i, -A->data[j][i]);\
        }\
    }\
    for(int i=0; i<a->nrow; i++){\
        for(int j=0; j<I->ncol; j++){\
            a->data[i][j] = A->data[i][a->ncol+b->ncol+j];\
        }\
        for(int j=0; j<b->ncol; j++){\
            b->data[i][j] = A->data[i][a->ncol+j];\
        }\
    }\
}\

NR__gaussj(float);



#undef SWAP
#undef NRANSI
