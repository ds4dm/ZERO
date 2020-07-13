/* #############################################
 *             This file is part of
 *                    ZERO
 *
 *             Copyright (c) 2020
 *     Released under the Creative Commons
 *        Zero v1.0 Universal License
 *
 *              Find out more at
 *        https://github.com/ds4dm/ZERO
 * #############################################*/


#include "support/utils.h"
#include <armadillo>
#include <boost/log/trivial.hpp>
#include <fstream>


arma::sp_mat
Utils::resizePatch(const arma::sp_mat &mat, const unsigned int nR, const unsigned int nC) {
  /**
 @brief Armadillo patch for resizing arma::sp_mat
 @details Armadillo arma::sp_mat::resize() is not robust as it initializes garbage
 values to new columns. This fixes the problem by creating new columns with
 guaranteed zero values. For arma::sp_mat
 */
  arma::sp_mat mMat(nR, nC);
  mMat.zeros();
  if (nR >= mat.n_rows && nC >= mat.n_cols) {
	 if (mat.n_rows >= 1 && mat.n_cols >= 1)
		mMat.submat(0, 0, mat.n_rows - 1, mat.n_cols - 1) = mat;
  } else {
	 if (nR <= mat.n_rows && nC <= mat.n_cols)
		mMat = mat.submat(0, 0, nR, nC);
	 else
		throw ZEROException(ZEROErrorCode::OutOfRange,
								  "Either both dimension should be smaller or larger.");
  }
  return mMat;
}

// For arma::mat
arma::mat Utils::resizePatch(const arma::mat &mat, const unsigned int nR, const unsigned int nC) {
  /**
 @brief Armadillo patch for resizing mat
 @details Armadillo mat::resize() is not robust as it initializes garbage
 values to new columns. This fixes the problem by creating new columns with
 guaranteed zero values. For arma::mat
 */
  arma::mat mMat(nR, nC);
  mMat.zeros();
  if (nR >= mat.n_rows && nC >= mat.n_cols) {
	 if (mat.n_rows >= 1 && mat.n_cols >= 1)
		mMat.submat(0, 0, mat.n_rows - 1, mat.n_cols - 1) = mat;
  } else {
	 if (nR <= mat.n_rows && nC <= mat.n_cols)
		mMat = mat.submat(0, 0, nR, nC);
	 else
		throw ZEROException(ZEROErrorCode::OutOfRange,
								  "Either both dimension should be smaller or larger.");
  }
  return mMat;
}

// For arma::vec
arma::vec Utils::resizePatch(const arma::vec &mat, const unsigned int nR) {
  /**
 @brief Armadillo patch for resizing arma::vec
 @details Armadillo arma::vec::resize() is not robust as it initializes garbage
 values to new columns. This fixes the problem by creating new columns with
 guaranteed zero values. For arma::vec
 */
  arma::vec mMat(nR);
  mMat.zeros();
  if (mat.n_rows > 0) {
	 if (nR > mat.n_rows)
		mMat.subvec(0, mat.n_rows - 1) = mat;
	 else
		mMat = mat.subvec(0, nR);
  }
  return mMat;
}

void Utils::appendSave(const arma::sp_mat &matrix, ///< The arma::sp_mat to be saved
							  const std::string   out,    ///< File name of the output file
							  const std::string   header, ///< A header that might be used to
																	///< check data correctness
							  bool erase                  ///< Should the matrix be appended to the
																	///< current file or overwritten
							  )
/**
 * Utility to append an arma::sp_mat to a data file.
 */
{
  // Using C++ file operations to copy the data into the target given by @out
  unsigned int nR{0}, nC{0}, nnz{0};

  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);

  nR  = matrix.n_rows;
  nC  = matrix.n_cols;
  nnz = matrix.n_nonzero;

  outfile << header << "\n";
  outfile << nR << "\t" << nC << "\t" << nnz << "\n";
  for (auto it = matrix.begin(); it != matrix.end(); ++it)
	 outfile << it.row() << "\t" << it.col() << "\t" << (*it)
				<< "\n"; // Write the required information of arma::sp_mat
  outfile << "\n";
  outfile.close(); // and close it
}

long int Utils::appendRead(arma::sp_mat &matrix, ///< Read and store the solution in this matrix.
									const std::string in, ///< File to read from (could be file very many
																 ///< data is appended one below another)
									long int pos, ///< Position in the long file where reading should start
									const std::string header ///< Any header to check data sanctity
									)
/**
 * Utility to read an arma::sp_mat from a long file.
 * @returns The end position from which the next data object can be read.
 */
{
  unsigned int nR = 0, nC = 0, nnz = 0;

  std::ifstream infile(in, std::ios::in);
  infile.seekg(pos);

  std::string headerCheckwith;
  infile >> headerCheckwith;

  if (header != "" && header != headerCheckwith)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Wrong header. Expected " + header + " found " + headerCheckwith);

  infile >> nR >> nC >> nnz;
  if (nR == 0 || nC == 0)
	 matrix.set_size(nR, nC);
  else {
	 arma::umat locations(2, nnz);
	 arma::vec  values(nnz);

	 unsigned int r = 0, c = 0;
	 double       val = 0;

	 for (unsigned int i = 0; i < nnz; ++i) {
		infile >> r >> c >> val;
		locations(0, i) = r;
		locations(1, i) = c;
		values(i)       = val;
	 }
	 matrix = arma::sp_mat(locations, values, nR, nC);
  }

  pos = infile.tellg();
  infile.close();

  return pos;
}

void appendSave(const std::vector<double> v,
					 const std::string         out,
					 const std::string         header,
					 bool                      erase) {
  /**
	* Utility to append an std::vector<double> to a data file.
	*/
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << header << "\n" << v.size() << "\n";
  for (const double x : v)
	 outfile << x << "\n";
  outfile.close();
}

long int
appendRead(std::vector<double> &v, const std::string in, long int pos, const std::string header) {
  unsigned long int size = 0;
  std::ifstream     infile(in, std::ios::in);
  infile.seekg(pos);
  /**
	* Utility to read an std::vector<double> from a long file.
	* @returns The end position from which the next data object can be read.
	*/

  std::string headerCheckwith;
  infile >> headerCheckwith;

  if (header != "" && header != headerCheckwith)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Wrong header. Expected " + header + " found " + headerCheckwith);

  infile >> size;

  v.resize(size);
  for (unsigned int i = 0; i < size; ++i)
	 infile >> v[i];
  pos = infile.tellg();
  infile.close();
  return pos;
}

void Utils::appendSave(const arma::vec & matrix, ///< The arma::vec to be saved
							  const std::string out,    ///< File name of the output file
							  const std::string header, ///< A header that might be used to
																 ///< check data correctness
							  bool erase                ///< Should the arma::vec be appended to the
																 ///< current file or overwritten
) {
  /**
	* Utility to append an arma::vec to a data file.
	*/
  // Using C++ file operations to copy the data into the target given by @out
  unsigned int nR{0};

  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);

  nR = matrix.n_rows;

  outfile << header << "\n";

  outfile << nR << "\n";
  for (auto it = matrix.begin(); it != matrix.end(); ++it)
	 outfile << (*it) << "\n"; // Write the required information of arma::sp_mat
  outfile << "\n";
  outfile.close(); // and close it
}

long int Utils::appendRead(arma::vec &matrix,    ///< Read and store the solution in this matrix.
									const std::string in, ///< File to read from (could be file very many
																 ///< data is appended one below another)
									long int pos, ///< Position in the long file where reading should start
									const std::string header ///< Any header to check data sanctity
) {
  /**
	* Utility to read an arma::vec from a long file.
	* @returns The end position from which the next data object can be read.
	*/
  unsigned int  nR;
  std::string   buffers;
  std::string   checkwith;
  std::ifstream inFile(in, std::ios::in);
  inFile.seekg(pos);

  inFile >> checkwith;
  if (header != "" && checkwith != header)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Wrong header. Expected " + header + " found " + checkwith);
  inFile >> nR;
  matrix.zeros(nR);
  for (unsigned int i = 0; i < nR; ++i) {
	 double val;
	 inFile >> val;
	 matrix.at(i) = val;
  }

  pos = inFile.tellg();
  inFile.close();

  return pos;
}

void Utils::appendSave(const long int    v,
							  const std::string out,
							  const std::string header,
							  bool              erase)
/**
 * Utility to save a long int to file
 */
{
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << header << "\n";
  outfile << v << "\n";
  outfile.close();
}

long int
Utils::appendRead(long int &v, const std::string in, long int pos, const std::string header) {
  /**
	* Utility to read a long int from a long file.
	* @returns The end position from which the next data object can be read.
	*/
  std::ifstream infile(in, std::ios::in);
  infile.seekg(pos);

  std::string headerCheckwith;
  infile >> headerCheckwith;

  if (header != "" && header != headerCheckwith)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Wrong header. Expected " + header + " found " + headerCheckwith);

  long int val;
  infile >> val;
  v = val;

  pos = infile.tellg();
  infile.close();

  return pos;
}

void Utils::appendSave(const unsigned int v,
							  const std::string  out,
							  const std::string  header,
							  bool               erase)
/**
 * Utility to save a long int to file
 */
{
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << header << "\n";
  outfile << v << "\n";
  outfile.close();
}

long int
Utils::appendRead(unsigned int &v, const std::string in, long int pos, const std::string header) {
  std::ifstream infile(in, std::ios::in);
  infile.seekg(pos);

  std::string headerCheckwith;
  infile >> headerCheckwith;

  if (header != "" && header != headerCheckwith)
	 throw ZEROException(ZEROErrorCode::InvalidData,
								"Wrong header. Expected " + header + " found " + headerCheckwith);

  unsigned int val;
  infile >> val;
  v = val;

  pos = infile.tellg();
  infile.close();

  return pos;
}

void Utils::appendSave(const std::string v, const std::string out, bool erase)
/**
 * Utility to save a long int to file
 */
{
  std::ofstream outfile(out, erase ? std::ios::out : std::ios::app);
  outfile << v << "\n";
  outfile.close();
}

long int Utils::appendRead(std::string &v, const std::string in, long int pos) {
  /**
	* Utility to read a std::string from a long file.
	* @returns The end position from which the next data object can be read.
	*/
  std::ifstream infile(in, std::ios::in);
  infile.seekg(pos);

  std::string val;
  infile >> val;
  v = val;

  pos = infile.tellg();
  infile.close();

  return pos;
}
unsigned long int Utils::vecToNum(std::vector<short int> binary) {
  unsigned long int number = 0;
  unsigned int      posn   = 1;
  while (!binary.empty()) {
	 short int bit = (binary.back() + 1) / 2; // The least significant bit
	 number += (bit * posn);
	 posn *= 2;         // Update place value
	 binary.pop_back(); // Remove that bit
  }
  return number;
}

std::vector<short int> Utils::numToVec(unsigned long int number, const unsigned long nCompl) {
  std::vector<short int> binary{};
  for (unsigned int vv = 0; vv < nCompl; vv++) {
	 binary.push_back(number % 2);
	 number /= 2;
  }
  std::for_each(binary.begin(), binary.end(), [](short int &vv) { vv = (vv == 0 ? -1 : 1); });
  std::reverse(binary.begin(), binary.end());
  return binary;
}

bool Utils::containsConstraint(const arma::sp_mat &A,
										 const arma::vec &   b,
										 const arma::vec &   lhs,
										 const double &      rhs,
										 const double        tol) {
  if (lhs.size() != A.n_cols)
	 return false;
  for (int i = 0; i < A.n_rows; ++i) {
	 bool res = true;
	 for (int j = 0; j < A.n_cols; ++j) {
		if (std::abs(lhs.at(j) - A.at(i, j)) > tol) {
		  res = false;
		  break;
		}
	 }
	 if (res && std::abs(b.at(i) - rhs) < tol) {
		return true;
	 }
  }
  return false;
}

bool Utils::containsElement(const arma::vec &b, const double &element, const double tol) {
  for (unsigned int i = 0; i < b.size(); ++i) {
	 if (std::abs(b.at(i) - element) < tol)
		return true;
  }
  return false;
}

bool Utils::containsRow(const arma::sp_mat &A, const arma::vec &row, const double tol) {

  if (row.size() != A.n_cols)
	 return false;
  for (int i = 0; i < A.n_rows; ++i) {
	 bool res = true;
	 for (int j = 0; j < A.n_cols; ++j) {
		if (std::abs(row.at(j) - A.at(i, j)) > tol) {
		  res = false;
		  break;
		}
	 }
	 if (res)
		return true;
  }
  return false;
}
bool Utils::containsConstraint(const arma::sp_mat &A,
										 const arma::vec &   b,
										 const arma::sp_mat &lhs,
										 const double &      rhs,
										 const double        tol) {
  if (lhs.n_rows > 1)
	 return false;
  arma::vec Ai = arma::vec{lhs};
  return Utils::containsConstraint(A, b, Ai, rhs, tol);
}


bool Utils::isZero(arma::mat M, double tol) noexcept {
  /**
	* @brief
	* Checking if a given matrix M is a zero matrix
	*
	* @p tol Tolerance, below which a number is treated as 0
	* @warning Tolerance < 0 always returns @p false with no error.
	*
	*/
  return ((abs(M).max() <= tol));
}

bool Utils::isZero(arma::sp_mat M, double tol) noexcept {
  /**
	* @brief
	* Checking if a given sparse matrix M is a zero matrix
	*
	* @p tol Tolerance, below which a number is treated as 0
	*
	*/
  if (M.n_nonzero != 0)
	 return false;

  return ((abs(M).max() <= tol));
}