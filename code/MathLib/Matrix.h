#pragma once

#pragma warning( disable : 4996)

#include <eigen3.2.5/Eigen/Dense>
#include <eigen3.2.5/Eigen/Sparse>
#include <cstdio>

typedef Eigen::Vector3d Vector3d;
typedef Eigen::VectorXd dVector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::Matrix2d Matrix2x2;
typedef Eigen::Matrix3d Matrix3x3;
typedef Eigen::Matrix4d Matrix4x4;
typedef Eigen::SparseMatrix<double> SparseMatrix;
typedef Eigen::Triplet<double> MTriplet;

inline void resize(SparseMatrix& sm, int rows, int cols) {
	if (sm.rows() != rows || sm.cols() != cols)
		sm.resize(rows, cols);
	sm.setZero();
}

inline void resize(dVector& v, int n) {
	if (v.size() != n)
		v.resize(n);
	v.setZero();
}

inline void resize(Matrix& m, int rows, int cols) {
	if (m.rows() != rows || m.cols() != cols)
		m.resize(rows, cols);
	m.setZero();
}

//TODO: want to write only upper (or lower?) values for symmetric matrices
template<class MATType>
void writeSparseMatrixDenseBlock(SparseMatrix& hes, int startX, int startY, const MATType& block, bool writeOnlyLowerDiagonalValues = false) {
    for (int i = 0; i < block.rows(); i++)
        for (int j = 0; j < block.cols(); j++)
            if (startX + i >= startY + j || !writeOnlyLowerDiagonalValues)
                hes.coeffRef(startX + i, startY + j) = block(i, j);
}

template<class MATType>
void writeSparseMatrixDenseBlockAdd(SparseMatrix& hes, int startX, int startY, const MATType& block, bool writeOnlyLowerDiagonalValues = false) {
    for (int i = 0; i < block.rows(); i++)
        for (int j = 0; j < block.cols(); j++)
            if (startX + i >= startY + j || !writeOnlyLowerDiagonalValues)
                hes.coeffRef(startX + i, startY + j) += block(i, j);
}

template<class MATType>
void addSparseMatrixDenseBlockToTriplet(std::vector<MTriplet>& triplets, int startX, int startY, const MATType& block, bool writeOnlyLowerDiagonalValues = false) {
    for (int i = 0; i < block.rows(); i++)
        for (int j = 0; j < block.cols(); j++)
        if (startX + i >= startY + j || !writeOnlyLowerDiagonalValues)
            triplets.push_back(MTriplet(startX + i, startY + j, block(i, j)));
}

//if the element at (row, col) is above the diagonal, it is skipped
inline void addMTripletToList_ignoreUpperElements(std::vector<MTriplet>& triplets, int row, int col, double val) {
	if (row >= col)
		triplets.push_back(MTriplet(row, col, val));
}

//if the element at (row, col) is above the diagonal, it is reflected on the lower diagonal
inline void addMTripletToList_reflectUpperElements(std::vector<MTriplet>& triplets, int row, int col, double val) {
	if (row >= col)
		triplets.push_back(MTriplet(row, col, val));
	else
		triplets.push_back(MTriplet(col, row, val));
}

//the element at (row, col) is mirrored, so it will be written symmetrically above and below the diagonal
inline void addMTripletToList_mirror(std::vector<MTriplet>& triplets, int row, int col, double val) {
	if (row == col) {
		triplets.push_back(MTriplet(row, col, val));
	}
	else {
		triplets.push_back(MTriplet(row, col, val));
		triplets.push_back(MTriplet(col, row, val));
	}
}

//write out the element as it comes
inline void addMTripletToList(std::vector<MTriplet>& triplets, int row, int col, double val) {
	triplets.push_back(MTriplet(row, col, val));
}

template <class MATType>
void print(char* fName, const MATType& mat) {
	FILE* fp; fp = fopen(fName, "w");

	for (int i = 0;i<mat.rows();i++) {
		for (int j = 0;j<mat.cols();j++) {
			fprintf(fp, "%f\t", mat.coeff(i, j));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}


#define ADD_HES_ELEMENT(list, i, j, v, w) addMTripletToList_reflectUpperElements(list, i, j, (v)*(w)) 
