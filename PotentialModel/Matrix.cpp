#include "stdafx.h"
#include <math.h>
#include "Matrix.h"


Matrix::Matrix(int n)
{
	m_n = N_PARAM;	
}

Matrix::~Matrix()
{
}

void Matrix::m_SetMatrix(double A[N_PARAM][N_PARAM])
{
	
	for (int i = 0; i < m_n; i ++) {
		for (int j = 0; j < m_n; j ++) {
			m_A[i][j] = A[i][j];
		}
	}
}

void Matrix::m_GetMatrix(double A[N_PARAM][N_PARAM])
{
	for (int i = 0; i < m_n; i ++) {
		for (int j = 0; j < m_n; j ++) {
			A[i][j] = m_A[i][j];
		}
	}
}

void Matrix::m_SetVector(double* b)
{
	for (int i = 0; i < m_n; i ++) {
		m_b[i] = b[i];
	}
}

void Matrix::m_GetVector(double* b)
{
	for (int i = 0; i < m_n; i ++) {
		b[i] = m_b[i];
	}
}

double Matrix::m_TriangDecomp()
{
// ���ò���ѡ��Ԫ��˹��ȥ����ʵ���� A �������Ƿֽ�
// ������� A �������� Cond(A) �Ľ���ֵ

	double ek, t, Anorm, Ynorm, Znorm, Cond;
	short i, j, k, l, m;

	double* work = new double [m_n];
	for (i = 0; i < m_n; i ++) {
		work[i] = 0.;
	}

	m_Pivot[m_n - 1] = 1;

	if (m_n == 1) {
		if (m_A[0][0] != 0.) {
			Cond = 1.;
		}
		else {
			Cond = 1.e+32;
		}
		delete []work;
		return Cond;
	}

//	���� A �����������������ֵ
	Anorm = 0.;
	for (j = 0; j < m_n; j ++) {
		t = 0.;
		for (i = 0; i < m_n; i ++) {
			t += fabs(m_A[i][j]);
		}
		if (t >= Anorm) Anorm = t;
	}

//	���ò���ѡ��Ԫ��˹��ȥ������ A ��Ϊ�����Ǿ���
	for (k = 0; k < m_n - 1; k ++) {
//		ѡ��Ԫ
		m = k;
		for (i = k + 1; i < m_n; i ++) {
			if (fabs(m_A[i][k]) > fabs(m_A[m][k])) m = i;
		}
		m_Pivot[k] = m;
		if (m != k) m_Pivot[m_n - 1] = -m_Pivot[m_n - 1];
//		����Ԫ�����Խ�����
		t = m_A[m][k];
		m_fswap(m_A[m][k], m_A[k][k]);
		if (t != 0.) {
//			���㱶������
			for (i = k + 1; i < m_n; i ++) {
				m_A[i][k] = -m_A[i][k] / t;
			}
//			����Ԫ������֮���������ʩ����ѡ��Ԫ��ȥ������ͬ������
			for (j = k + 1; j < m_n; j ++) {
				t = m_A[m][j];
				m_fswap(m_A[m][j], m_A[k][j]);
				if (t == 0.) continue;
				for (i = k + 1; i < m_n; i ++) {
					m_A[i][j] += m_A[i][k] * t;
				}
			}
		}
	}

//	��ⷽ���� A'*y=e������ A' Ϊ A ��ת�þ���e �Ƿ���Ϊ��1 ������
	for (k = 0; k < m_n; k ++) {
		t = 0.;
		if (k > 0) {
			for (i = 0; i < k; i ++) {
				t += m_A[i][k] * work[i];
			}
		}
		if (t < 0.) {
			ek = -1.;
		}
		else {
			ek = 1.;
		}
		if (m_A[k][k] == 0.) {
			Cond = 1.E+32;
			delete []work;
			return Cond;
		}
		work[k] = -(ek + t) / m_A[k][k];
	}

	for (l = 0; l < m_n - 1; l ++) {
		k = m_n - l - 2;
		t = 0.;
		for (i = k + 1; i < m_n; i ++) {
			t += m_A[i][k] * work[k];
		}
		work[k] = t;
		m = m_Pivot[k];
		if (m != k)	{
			m_fswap(work[m], work[k]);
		}
	}

//	���� ��y��
	Ynorm = 0.;
	for (i = 0; i < m_n; i ++) {
		Ynorm += fabs(work[i]);
	}

//	��ⷽ���� A*z=y
	m_SetVector(work);
	m_BackSubstitute();
	m_GetVector(work);

//	���� ��z��
	Znorm = 0.;
	for (i = 0; i < m_n; i ++) {
		Znorm += fabs(work[i]);
	}

	
//	���� Cond(A) �Ľ���ֵ
	Cond = Anorm * Znorm / Ynorm;
	if (Cond < 1.) Cond = 1.;
	delete []work;
	return Cond;
}

void Matrix::m_BackSubstitute()
{
//	������Է����� A*X=b������ A �Ѿ���˹��ȥ����Ϊ�����Ǿ���

	int i, k, l, m;
	double t;

//	���Ҷ��� b ʩ����ѡ��Ԫ��ȥ������ͬ������
	if (m_n > 1) {
		for (k = 0; k < m_n - 1; k ++) {
			m = m_Pivot[k];
			t = m_b[m];
			m_fswap(m_b[m], m_b[k]);
			for (i = k + 1; i < m_n; i ++) {
				m_b[i] += m_A[i][k] * t;
			}
		}

//		�ش�����÷�����Ľ�
		for (l = 0; l < m_n - 1; l ++) {
			k = m_n - l - 1;
			m_b[k] /= m_A[k][k];
			t = -m_b[k];
			for (i = 0; i < m_n - l - 1; i ++) {
				m_b[i] += m_A[i][k] * t;
			}
		}
	}
	m_b[0] /= m_A[0][0];
}

void Matrix::m_fswap(double& a, double& b)
{
	double t;
	t = a;
	a = b;
	b = t;
}
