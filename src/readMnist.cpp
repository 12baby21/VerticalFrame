#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

// 改变字节序
int reverseEndian(int i)
{
	unsigned char ch1, ch2, ch3, ch4;
	ch1 = i & 0xFF;
	ch2 = (i >> 8) & 0xFF;
	ch3 = (i >> 16) & 0xFF;
	ch4 = (i >> 24) & 0xFF;
	return ((int)ch1 << 24) + ((int)ch2 << 16) + ((int)ch3 << 8) + ch4;
}

bool readMnistLabel(const string filename, vector<char> &labels)
{
	ifstream file(filename, ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int num_images = 0;
		file.read((char *)&magic_number, sizeof(magic_number));
		file.read((char *)&num_images, sizeof(num_images));
		magic_number = reverseEndian(magic_number);
		num_images = reverseEndian(num_images);

		for (int i = 0; i < num_images; i++)
		{
			char label = 0;
			file.read(&label, sizeof(char));
			labels.push_back(label);
		}
	}
	else
	{
		cout << "文件打开失败" << endl;
		return false;
	}
	return true;
}

bool readMnistImages(const string filename, vector<vector<float>> &images)
{
	ifstream file(filename, ios::binary);
	if (file.is_open())
	{
		int magic_number = 0;
		int num_images = 0;
		int row = 0;
		int column = 0;
		unsigned char label;
		file.read((char *)&magic_number, sizeof(magic_number));
		file.read((char *)&num_images, sizeof(num_images));
		file.read((char *)&row, sizeof(row));
		file.read((char *)&column, sizeof(column));
		magic_number = reverseEndian(magic_number);
		num_images = reverseEndian(num_images);
		row = reverseEndian(row);
		column = reverseEndian(column);

		for (int i = 0; i < num_images; ++i)
		{
			vector<float> tmp;
			for (int j = 0; j < row * column; ++j)
			{
				unsigned char image = 0;
				file.read((char *)&image, sizeof(char));
				tmp.push_back(image);
			}
			images.push_back(tmp);
		}
	}
	else
	{
		cout << "文件打开失败" << endl;
		return false;
	}
	return true;
}

void processMnistLabeltoOddandEven(vector<vector<long>> &newLabels, vector<char> &labels)
{
	int n = labels.size();
	for (int i = 0; i < n; ++i)
	{
		labels[i] = labels[i] % 2;
		if (labels[i] == 1)
		{
			newLabels[i][1] = 1;
		}
		else
		{
			newLabels[i][0] = 1;
		}
	}
}

void normalizeMnistImage(vector<vector<float>> &images)
{
	int n = images.size();
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < 784; ++j)
		{
			images[i][j] = images[i][j] / 255.0;
		}
	}
}

/**
 * @brief 对数据集进行特征分割
 *
 * @param A 参与方A所拥有的特征
 * @param B 参与方B所拥有的特征
 * @param images 原始特征
 */
void partitionFeatures(vector<vector<float>> &A, vector<vector<float>> &B, vector<vector<float>> &images)
{
	int num_image = images.size();
	int featureNum = images[0].size();
	for (int i = 0; i < num_image; ++i)
	{
		for (int j = 0; j < featureNum / 2; ++j)
		{
			A[i].push_back(images[i][j]);
			B[i].push_back(images[i][featureNum + j]);
		}
	}
}