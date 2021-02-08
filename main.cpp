#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>
#include <random>
#include <climits>
#include <cmath>
#include <chrono>
#include <conio.h>

#define PERFOMANCE 3732764

using namespace std;

// Функция, проверяющая является ли n степенью k.
bool checkPow(int n, int k) {
	if(n < k)
		return false;
	while(n >= 1) {
		if(n == 1)
			return true;
		n/=k;
	}
	return false;
}
// Структура узла, используемая классом Tree.
struct Node {
	unsigned long L, M, item;
	Node *parent;
	vector<Node *> sons;
	// Конструктор листа
	Node(unsigned long i, Node *pnt = nullptr): item(i), parent(pnt) {}
	// Конструктор внутреннего узла
	Node(Node *l_son, Node *m_son, Node *pnt = nullptr) {
		L = MAX(l_son);
		M = MAX(m_son);
		parent = pnt;
		sons.push_back(l_son);
		sons.push_back(m_son);
	}
	// Функция проверки, является ли вызывающий узел листом.
	bool isLeaf() {
		return sons.empty();
	}
	// Функция собирает элементы дерева, корнем которого является узел node,
	// в вектор и возвращает его на выходе.
	static vector<unsigned long> toVector(Node *node) {
		vector<unsigned long> v;
		if(node->isLeaf()) {
			v.push_back(node->item);
		}
		else {
			vector<unsigned long> temp;
			for(int i = 0; i < node->sons.size(); i++) {
				temp = toVector(node->sons[i]);
				v.insert(v.end(), temp.begin(), temp.end());
			}
		}
		return v;
	}
	// Функция возвращает строку, содержащую структуру дерева,
	// корнем которого является вызывающий узел. Параметр height
	// определяет полную высоту дерева.
	string toString(int height, string spaces = "\t\t\t\t") {
		string res = "";
		for(int i = this->HEIGHT(); i < height; i++)
			res += spaces;
		res += "    " + this->nodeToString() + '\n';
		for(int i = 0; i < 100; i++)
			res += '-';
		res += '\n';
		for(int i = 0; i < this->sons.size(); i++)
			res += this->sons[i]->toString(height);
		return res;
	}
	// Функция проверки вхождения элемента item в вызывающем узле
	// или его сыновьях.
	bool isContained(unsigned long item) {
		if(this->isLeaf())
			return item == this->item;

		bool res = false;
		for(int i = 0; i < sons.size(); i++)
			res = res || sons[i]->isContained(item);

		return res;
	}
	// Функция преобразования узла в информативную строку.
	// Для листа возвращается знаение, хранящееся в нем, а для
	// внутреннего узла его значения L и M.
	string nodeToString() {
		return this->isLeaf()?to_string(this->item):
			(to_string(this->L) + '-' + to_string(this->M));
	}
	// Функция поиска максимального элемента в поддереве с корнем node.
	unsigned long MAX(Node *node) {
		Node *u = node;
		for(; !u->isLeaf(); u = u->sons[u->sons.size() - 1]);
		return u->item;
	}
	// Функция поиска минимального элемента в поддереве с корнем node.
	unsigned long MIN(Node *node) {
		Node *u = node;
		for(; !u->isLeaf(); u = u->sons[0]);
		return u->item;
	}
	// Функция возвращает высоту вызывающего узла.
	int HEIGHT() {
		int height;
		Node *u;
		for(u = this, height = 0; !u->isLeaf(); u = u->sons[0], height++);
		return height;
	}
};
// Класс, реализующий 2-3 дерево.
class Tree {
	// Функция поиска подходящего узла для вставки элемента item.
	// Используется исключительно самим объектом класса.
	Node *getNodeForItem(unsigned long item) {
		Node *current = root;
		if(!current->isLeaf())
			// Если первый сын не является листом, то
			// не являются и остальные
			while(!current->sons[0]->isLeaf()) {
				if(item <= current->L)
					current = current->sons[0];
				else if(item <= current->M)
					current = current->sons[1];
				else if(current->sons.size() == 2){
					current->M = item;
					current = current->sons[1];
				}
				else current = current->sons[2];
			}
		return current;
	}
public:
	Node *root;
	//  Конструктор пустого дерева.
	Tree() {
		root = nullptr;
	}
	// Конструктор дерева, полученного в результате объединения
	// деревьев T1 и T2.
	Tree(Tree *T1, Tree *T2) {
		root = t_combine(T1, T2)->root;
	}
	// Конструктор дерева из элементов вектора v.
	Tree(vector<unsigned long> v) {
		root = nullptr;
		t_addItems(v);
	}
	// Конструктор дерева критической организации по заданным параметрам.
	// n - количество элементов, start_value - начальное(минимальное) значение,
	// tree_type флаг вида критической организации (2 - каждый узел имеет по 2
	// сына; 3 - каждый узел имеет по 3 сына, n должен быть степенью 3!!!). Деревья
	// заполняются n значениями, каждое из которых отличается друг от друга на 2*k,
	// где k принадлежит [0; n-1]. Т.е. если начальное значение четное - дерево
	// заполнится четными числами, иначе нечетными.
	Tree(int n, unsigned long start_value, int tree_type) {
		root = nullptr;
		vector<unsigned long> v;
		unsigned long value = start_value;

		if(tree_type == 2) {
			for(int i = 0; i < n; i++, value+=2)
				v.push_back(value);

			t_addItems(v);
		}
		else if(tree_type == 3) {
			if(checkPow(n, 3)) {
				if(n == 3) {
					v.push_back(value);
					v.push_back(value + 2);
					v.push_back(value + 4);

					t_addItems(v);
				}
				else {
					int part_size = n/3;
					Tree *T1, *T2, *T3;

					T1 = new Tree(part_size, start_value, tree_type);
					T2 = new Tree(part_size, start_value + 2 * part_size, tree_type);
					T3 = new Tree(part_size, start_value + 4 * part_size, tree_type);

					root = combine(combine(T1, T2), T3)->root;
				}
			}
		}
	}
	// Функция поиска максимального элемента в вызывающем дереве.
	unsigned long MAX() {
		Node *u = root;
		for(; !u->isLeaf(); u = u->sons[u->sons.size() - 1]);
		return u->item;
	}
	// Функция поиска минимального элемента в вызывающем дереве.
	unsigned long MIN() {
		Node *u = root;
		for(; !u->isLeaf(); u = u->sons[0]);
		return u->item;
	}
	// Функция возвращает высоту вызывающего дерева.
	int HEIGHT() {
		if(root == nullptr)
			return 0;
		int height;
		Node *u;
		for(u = root, height = 0; !u->isLeaf(); u = u->sons[0], height++);
		return height;
	}
	// Функция выводит вызывающее дерево.
	void showTree() {
		cout << endl;
		for(int i = 0; i < 100; i++)
			cout << '=';
		cout << endl << this->root->toString(this->HEIGHT());
		for(int i = 0; i < 100; i++)
			cout << '=';
		cout << endl;
	}
	// Функция собирает элементы вызывающего дерева в вектор и возвращает
	// его на выходе.
	vector<unsigned long> toVector() {
		return Node::toVector(root);
	}
	// Процедура балансировки узла node. Используется в случаях, когда
	// узел node имеет 4 сыновей.
	void balance(Node *node) {
		// Если node не имеет родителя - создать его. После чего сделать
		// его двух правых сыновей сыновьями его ближайшего правого брата (брат создается).
		if(node->parent == nullptr) {
			Node *sec_node = new Node(node->sons[2], node->sons[3]);
			node->sons[2]->parent = sec_node;
			node->sons[3]->parent = sec_node;
			node->sons.erase(node->sons.end() - 2, node->sons.end());

			Node *new_root = new Node(node, sec_node);
			node->parent = new_root;
			sec_node->parent = new_root;
			root = new_root;
		}
		// Иначе сразу сделать двух правых сыновей узла node сыновьями
		// его ближайшего правого брата (брат создается).
		else {
			Node *sec_node = new Node(node->sons[2], node->sons[3], node->parent);
			node->sons[2]->parent = sec_node;
			node->sons[3]->parent = sec_node;
			node->sons.erase(node->sons.end() - 2, node->sons.end());
			for(int i = 0; i < node->parent->sons.size(); i++) {
				if(node->parent->sons[i] == node) {
					node->parent->sons.insert(node->parent->sons.begin() + i + 1, sec_node);
					// Корректировка значений L и M
					if(i == 0) {
						node->parent->L = root->MAX(node);
						node->parent->M = root->MAX(sec_node);
					}
					if(i == 1)
						node->parent->M = root->MAX(node);
					// Рекурсивная балансировка родительского узла
					if(node->parent->sons.size() == 4)
						balance(node->parent);
					return;
				}
			}
		}
	}
	// Процедура добавления элемента item в дерево
	void addItem(unsigned long item) {
		// Если дерево пустое, создаем его
		if (root == nullptr) {
			root = new Node(item);
			return;
		}
		// Находим подходящий для вставки узел
		Node *s_node = getNodeForItem(item);
		// Если элемент item уже содержится в дереве, прекращаем операцию
		if(s_node->isContained(item))
			return;
		// Создание нового узла-листа
		Node *new_node = new Node(item, s_node);
		// Если подходящий узел является листом (дерево состоит из 1 узла),
		// создаем новый корень и делаем имеющиеся узлы-листы его сыновьями
		if(s_node->isLeaf()) {
			Node *new_root;
			if(s_node->item < new_node->item)
				new_root = new Node(s_node, new_node);
			else
				new_root = new Node(new_node, s_node);
			root = new_root;
			new_node->parent = new_root;
			s_node->parent = new_root;
			return;
		}
		// Поиск места для вставки в подходящем узле и последующая вставка
		for(int i = 0; i < s_node->sons.size(); i++)
			if(new_node->item < s_node->sons[i]->item) {
				s_node->sons.insert(s_node->sons.begin() + i, new_node);
				// Корректировка значений L и M
				if(i == 0) {
					s_node->L = root->MAX(s_node->sons[0]);
					s_node->M = root->MAX(s_node->sons[1]);
				}
				if(i == 1)
					s_node->M = root->MAX(s_node->sons[1]);

				// Выполнение балансировки
				if(s_node->sons.size() == 4)
					balance(s_node);
				return;
			}
		s_node->sons.push_back(new_node);
		// Выполнение балансировки
		if(s_node->sons.size() == 4)
			balance(s_node);
	}
	// Процедура последовательной вставки элементов вектора v
	// в вызывающее дерево в выводом шкалы прогресса.
	void t_addItems(vector<unsigned long> v) {
		int bar_width = 100, progress = 0;
		if(v.size() < 10 * bar_width) {
			for(int i = 0; i < v.size(); i++)
				addItem(v[i]);
			return;
		}
		int segment_size = v.size() / bar_width;
		for(int i = 0; i < v.size(); i++){
			if(i % segment_size == 0) {
				cout.flush();
				cout << "|";
				for(int j = 0; j < bar_width; j++){
					if(j < progress)
						cout << "#";
					else
						cout << " ";
				}
				cout << "|" << int(((float)progress/bar_width) * 100) << "%\r";
				progress++;
			}
			addItem(v[i]);
		}
	}
	// Процедура последовательной вставки элементов вектора v
	// в вызывающее дерево.
	void addItems(vector<unsigned long> v) {
		for(int i = 0; i < v.size(); i++)
			addItem(v[i]);
	}
	// Функция combine(...) с выводом времени выполнения программы.
	static Tree *t_combine(Tree *T1, Tree *T2) {
		long double start_time, end_time;
		start_time = clock();
		Tree *T3 = combine(T1, T2);
		end_time = clock();
		cout << endl << "Runtime = " << (end_time - start_time)/1000.0 << " sec" << endl;
		return T3;
	}
	// Функция объединения деревьев, возвращает новое дерево.
	static Tree *combine(Tree *T1, Tree *T2) {
		unsigned long min1 = T1->MIN(), max1 = T1->MAX(),
							min2 = T2->MIN(), max2 = T2->MAX(),
							height1 = T1->HEIGHT(), height2 = T2->HEIGHT();
		Tree *new_Tree = (height1 < height2)?T2:T1;
		// Анализ диапазонов значений деревьев.
		if((min1 <= min2 && min2 <= max1) || (min1 <= max2 && max2 <= max1)) {
			// cout << endl << "Ranges overlap!" << endl;
			// Последовательная вставка элементов.
			vector<unsigned long> v;
			if(height1 < height2)
				v = T1->toVector();
			else 
				v = T2->toVector();
			new_Tree->addItems(v);
		}
		else {
			// cout << endl << "Ranges do not overlap!" << endl;
			// Объединение деревьев созданием общего узла.
			Node *u = new_Tree->root;
			// Если дерево с большим диапазоном находится по указателю T2.
			if(max1 < max2) {
				// Если высота дерева слева меньше высоты дерева справа.
				if(height1 < height2) {
					for(int i = 1; i < (height2 - height1); i++, u = u->sons[0]);
					u->sons.insert(u->sons.begin(), T1->root);
					T1->root->parent = u;
					u->M = u->L;
					u->L = max1;
					if(u->sons.size() == 4)
						new_Tree->balance(u);
				}
				// Иначе если высота дерева слева больше высоты дерева справа.
				else if(height1 > height2) {
					for(int i = 1; i < (height1 - height2); i++, u = u->sons[u->sons.size() - 1]);
					u->sons.push_back(T2->root);
					T2->root->parent = u;
					if(u->sons.size() == 4)
						new_Tree->balance(u);
				}
				// Иначе если высоты равны.
				else {
					new_Tree = new Tree();
					new_Tree->root = new Node(T1->root, T2->root);
					T1->root->parent = new_Tree->root;
					T2->root->parent = new_Tree->root;
				}
			}
			// Иначе если дерево с большим диапазоном находится по указателю T1.
			else {
				// Если высота дерева слева меньше высоты дерева справа.
				if(height1 < height2) {
					for(int i = 1; i < (height2 - height1); i++, u = u->sons[u->sons.size() - 1]);
					u->sons.push_back(T1->root);
					T1->root->parent = u;
					if(u->sons.size() == 4)
						new_Tree->balance(u);
				}
				// Иначе если высота дерева слева больше высоты дерева справа.
				else if(height1 > height2) {
					for(int i = 1; i < (height1 - height2); i++, u = u->sons[0]);
					u->sons.insert(u->sons.begin(), T2->root);
					T2->root->parent = u;
					u->M = u->L;
					u->L = max2;
					if(u->sons.size() == 4)
						new_Tree->balance(u);
				}
				// Иначе если высоты равны.
				else {
					new_Tree = new Tree();
					new_Tree->root = new Node(T2->root, T1->root);
					T1->root->parent = new_Tree->root;
					T2->root->parent = new_Tree->root;
				}
			}
		}
		return new_Tree;
	}
};
// Процедура тестирования: тест №1.
void test_1() {
	string tabs = "\t";
	cout << endl << tabs + "TEST 1:" << endl << endl;

	Tree *T1 = new Tree();
	T1->addItem(1), T1->addItem(2);

	Tree *T2 = new Tree();
	T2->addItem(3), T2->addItem(4);
	T2->addItem(5), T2->addItem(6);

	cout << tabs + "TREE 1:" << endl;
	T1->showTree();
	cout << endl << endl;
	getch();

	cout << tabs + "TREE 2:" << endl;
	T2->showTree();
	cout << endl << endl;
	getch();

	Tree *T3 = Tree::t_combine(T1, T2);

	cout << tabs + "RESULT:" << endl;
	T3->showTree();
	cout << endl << endl;
}
// Процедура тестирования: тест №2.
void test_2() {
	string tabs = "\t";
	cout << endl << tabs + "TEST 2:" << endl << endl;

	Tree *T1 = new Tree();
	T1->addItem(7), T1->addItem(9), T1->addItem(14);

	Tree *T2 = new Tree();
	T2->addItem(4), T2->addItem(6);

	cout << tabs + "TREE 1:" << endl;
	T1->showTree();
	cout << endl << endl;
	getch();

	cout << tabs + "TREE 2:" << endl;
	T2->showTree();
	cout << endl << endl;
	getch();

	Tree *T3 = Tree::t_combine(T1, T2);

	cout << tabs + "RESULT:" << endl;
	T3->showTree();
	cout << endl << endl;
}
// Процедура тестирования: тест №3.
void test_3() {
	string tabs = "\t";
	cout << endl << tabs + "TEST 3:" << endl << endl;

	Tree *T1 = new Tree();
	T1->addItem(3), T1->addItem(5), T1->addItem(7);

	Tree *T2 = new Tree();
	T2->addItem(9);

	cout << tabs + "TREE 1:" << endl;
	T1->showTree();
	cout << endl << endl;
	getch();

	cout << tabs + "TREE 2:" << endl;
	T2->showTree();
	cout << endl << endl;
	getch();

	Tree *T3 = Tree::t_combine(T1, T2);

	cout << tabs + "RESULT:" << endl;
	T3->showTree();
	cout << endl << endl;
}
// Процедура тестирования: тест №4.
void test_4() {
	string tabs = "\t";
	cout << endl << tabs + "TEST 4:" << endl << endl;

	Tree *T1 = new Tree();
	T1->addItem(5), T1->addItem(9);

	Tree *T2 = new Tree();
	T2->addItem(3), T2->addItem(4);
	T2->addItem(6), T2->addItem(8);

	cout << tabs + "TREE 1:" << endl;
	T1->showTree();
	cout << endl << endl;
	getch();

	cout << tabs + "TREE 2:" << endl;
	T2->showTree();
	cout << endl << endl;
	getch();

	Tree *T3 = Tree::t_combine(T1, T2);

	cout << tabs + "RESULT:" << endl;
	T3->showTree();
	cout << endl << endl;
}
// Функция измерения производительности компьютера. Возвращает
// количество выполненных операций инкремента за 3 секунды.
unsigned long measurePerfomance() {
	unsigned long i = 0;
	chrono::high_resolution_clock::time_point start, end;

	start = end = chrono::high_resolution_clock::now();
	for(; chrono::duration_cast<chrono::duration<double>>(end-start).count() < 1.0; 
		end = chrono::high_resolution_clock::now(), ++i);

	return i;
}
// Функция прогноза времени выполнения операции объединения деревьев T1 и T2.
vector<long double> predictTime(Tree *T1, Tree *T2) {
	int n, m;
	vector<long double> v;
	
	// Определение коэффициентов для выражений.
	long double min_C1 = 0.00002;
	long double avrg_C1 = 4.186062301413196e-07, avrg_C2 = 0.0003879197716936213,
		avrg_C3 = -3.034775198489289e-07, avrg_C4 = -4.993704794171108e-05;
	// long double avrg_C1 = 4.955516592959107e-07, avrg_C2 = 0.00047020803651980363,
	// 	avrg_C3 = 3.069586544567334e-06, avrg_C4 = -0.002014022810620488;
		
	long double max_C1 = 7.199413644755995e-07, max_C2 = 0.000164229897998686,
		max_C3 = 6.886464721002315e-06, max_C4 = -6.975900654599681e-05;
	// long double max_C1 = 5.023048361019243e-07, max_C2 = 0.00013684868387917233,
	// 	max_C3 = 6.643799855089625e-06, max_C4 = 2.9634135544018295e-05;
	// long double max_C1 = 5.027276855092602e-07, max_C2 = 0.0001682928073547816,
	// 	max_C3 = 6.743185074712246e-06, max_C4 = -0.00034931927744493054;

	if(T1->HEIGHT() >= T2->HEIGHT()) {
		n = T1->toVector().size();
		m = T2->toVector().size();
	}
	else {
		n = T2->toVector().size();
		m = T1->toVector().size();
	}

	if(n < 150 || m < 150)
		return v;

	// Получение коэффициента производительности исполняющей машины.
	unsigned long current_perfomance = measurePerfomance();

	cout << endl << current_perfomance << ' ' << PERFOMANCE << endl;

	// Вычисление и сохранение минимального, среднего и максимального времени
	// выполнения с учетом различия производительности машин.
	v.push_back((long double)(min_C1 * ((long double)PERFOMANCE / current_perfomance)));

	v.push_back((long double)(avrg_C1 * ((long double)log(n) / log(2.8)) * m + avrg_C2 * ((long double)log(n) / log(2.8)) + 
		avrg_C3 * m + avrg_C4) * ((long double)PERFOMANCE / current_perfomance));
	v.push_back((long double)(max_C1 * ((long double)log(n) / log(2)) * m + max_C2 * ((long double)log(n) / log(2)) + 
		max_C3 * m + max_C4) * ((long double)PERFOMANCE / current_perfomance));

	cout << "\n\nBest time = " << v[0] << endl;
	cout << "Average time = " << v[1] << endl;
	cout << "Worst time = " << v[2] << endl << endl;

	return v;
}
// Процедура тестирования с заданием параметров.
void newTest() {
	char cm;
	int size;
	double temp;
	vector<unsigned long> v;
	Tree *T1, *T2, *T3;
	chrono::high_resolution_clock::time_point start_time, end_time;

	cout << endl << "\tSETTING TREE 1:";

	cout << endl << "Press 1 to read the data for the tree from the file \"tree1_data.txt\"";
	cout << endl << "Press 2 to generate random tree.";
	cout << endl << "To exit press any other key!" << endl;

	cm = getch();
	if(cm == '1') {
		ifstream fin("tree1_data.txt");
		while(!fin.eof()) {
			fin >> temp;
			v.push_back(temp);
		}
		fin.close();

		T1 = new Tree(v);
	}
	else if(cm == '2') {
		cout << "Enter the number of items: ";
		cin >> size;

		mt19937 mersenne(time(NULL));

		for(int i = 0; i < size; i++)
			v.push_back(mersenne());

		T1 = new Tree(v);
	}
	else {
		cout << "Incorrect key. To continue press ENTER!" << endl;
		return;
	}

	v.clear();

	cout << endl << "Show the tree?(y/n)" << endl;
	cm = getch();
	if(cm == 'y')
		T1->showTree();

	cout << endl << "\tSETTING TREE 2:";
	
	cout << endl << "Press 1 to read the data for the tree from the file \"tree2_data.txt\"";
	cout << endl << "Press 2 to generate random tree.";
	cout << endl << "To exit press any other key!" << endl;

	cm = getch();
	if(cm == '1') {
		ifstream fin("tree2_data.txt");
		while(!fin.eof()) {
			fin >> temp;
			v.push_back(temp);
		}
		fin.close();

		T2 = new Tree(v);
	}
	else if(cm == '2') {
		cout << "Enter the number of items: ";
		cin >> size;

		mt19937 mersenne(time(NULL));

		for(int i = 0; i < size; i++)
			v.push_back(mersenne());

		T2 = new Tree(v);
	}
	else {
		cout << "Incorrect key. To continue press ENTER!" << endl;
		return;
	}

	cout << endl << "Show the tree?(y/n)" << endl;
	cm = getch();
	if(cm == 'y')
		T2->showTree();

	cout << endl << "\tECECUTION..." << endl;

	predictTime(T1, T2);

	start_time = chrono::high_resolution_clock::now();
	T3 = Tree::combine(T1, T2);
	end_time = chrono::high_resolution_clock::now();
	cout << "Runtime = ";
	cout << chrono::duration_cast<chrono::duration<long double>>(end_time - start_time).count();
	cout << endl;

	vector<unsigned long> nv = T3->toVector();
	cout << "The new tree contains " << nv.size() << " items" << endl;

	cout << endl << "Print the tree to console? (y/n)\n";

	cm = getch();
	if(cm == 'y')
		T3->showTree();

	cout << endl << "Write the tree to the file \"output.txt\"? (y/n) ";

	cm = getch();
	if(cm == 'y'){
		cout << "\n\nPlease wait...\t";
		string buff;
		buff = T3->root->toString(T3->HEIGHT());
		ofstream fout("output.txt");
		fout << buff;
		fout.close();
		cout << "Ready." << endl;
	}
	cout << "\n\nTo continue press ENTER!" << endl;

	return;
}
// Функция замера времени выполнения процедуры Tree::combine(...).
// Объединяемые деревья задаются в самой функции. Параметр n -
// количество выполнений каждого из тестов.
void timebest(int n = 10) {
	chrono::high_resolution_clock::time_point start_time, end_time;
	double average_time;
	vector<chrono::duration<double>> res;
	Tree *T1, *T2, *T3;
	int size1 = 16, size2 = 81;

	cout << "\n\tBest case testing.\n\n";

	cout << "\tTEST 1.\nn = " << size1 << ", m = " << size2 << ".\n";

	mt19937 mersenne(time(NULL));
	for(int i = 0; i < n; i++) {
		unsigned long m = mersenne();

		T1 = new Tree(size1, m, 2);
		T2 = new Tree(size2, m + size1 * 2, 3);

		start_time = chrono::high_resolution_clock::now();
		T3 = new Tree(T1, T2);
		end_time = chrono::high_resolution_clock::now();

		res.push_back(chrono::duration_cast<chrono::duration<double>>(end_time - start_time));

		// cout << "Round " << i + 1 << "\\" << n << ": ";
		// cout << res[res.size()].count() << " sec\n";
	}
	average_time = 0;
	for(int i = 0; i < res.size(); i++) {
		average_time += res[i].count();
	}
	average_time /= res.size();
	cout << "Average time = " << average_time << " sec\n\n";

	res.clear();
}
// Функция замера времени выполнения процедуры Tree::combine(...).
// Объединяемые деревья задаются в самой функции. Параметр n -
// количество выполнений каждого из тестов.
void timeworst(int n = 10) {
	chrono::high_resolution_clock::time_point start_time, end_time;
	double average_time;
	vector<chrono::duration<double>> res;
	Tree *T1, *T2, *T3;
	int size1 = 59049, size2 = 1024;

	cout << "\n\tWorst case testing.\n\n";

	cout << "\tTEST 1.\nn = " << size1 << ", m = " << size2 << ".\n";

	mt19937 mersenne(time(NULL));
	for(int i = 0; i < n; i++) {

		T1 = new Tree(size1, 1, 3);
		T2 = new Tree(size2, 2, 2);

		start_time = chrono::high_resolution_clock::now();
		T3 = new Tree(T1, T2);
		end_time = chrono::high_resolution_clock::now();

		res.push_back(chrono::duration_cast<chrono::duration<double>>(end_time - start_time));

		// cout << "Round " << i + 1 << "\\" << n << ": ";
		// cout << res[res.size()].count() << " sec\n";
	}
	average_time = 0;
	for(int i = 0; i < res.size(); i++) {
		average_time += res[i].count();
	}
	average_time /= res.size();
	cout << "Average time = " << average_time << " sec\n\n";

	res.clear();
}
// Функция замера времени выполнения процедуры Tree::combine(...).
// Объединяемые деревья задаются в самой функции. Параметр n -
// количество выполнений каждого из тестов.
void timerand(int n = 10) {
	chrono::high_resolution_clock::time_point start_time, end_time;
	double average_time;
	vector<chrono::duration<double>> res;
	Tree *T1, *T2, *T3;
	int size1 = 59049, size2 = 59049;
	vector<unsigned long> v1, v2;

	cout << "\n\tRand case testing.\n\n";

	cout << "\tTEST 1.\nn = " << size1 << ", m = " << size2 << ".\n";

	mt19937 mersenne(time(NULL));
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < size1; j++)
			v1.push_back(mersenne());
		for(int j = 0; j < size2; j++)
			v2.push_back(mersenne());

		T1 = new Tree(v1);
		T2 = new Tree(v2);

		v1.clear(), v2.clear();

		start_time = chrono::high_resolution_clock::now();
		T3 = new Tree(T1, T2);
		end_time = chrono::high_resolution_clock::now();

		res.push_back(chrono::duration_cast<chrono::duration<double>>(end_time - start_time));

		// cout << "Round " << i + 1 << "\\" << n << ": ";
		// cout << res[res.size()].count() << " sec\n";
	}
	// average_time = 0;
	for(int i = 0; i < res.size(); i++) {
		average_time += res[i].count();
	}
	average_time /= res.size();
	cout << "Average time = " << average_time << " sec\n\n";

	res.clear();
}
// Главная функция, реализующая меню выбора теста.
int main() {
	// timerand(50);
	// timeworst(50);
	// timebest(50);

	char cm;
	cout << endl << "\tMAIN MENU"; 
	cout << endl << "To start standard test enter its number(1-4).";
	cout << endl << "To create new test press 5";
	cout << endl << "To exit press \'q\'!" << endl;
	while((cm = getch()) != 'q') {
		switch (cm) {
			case '1':
				test_1();
				break;
			case '2':
				test_2();
				break;
			case '3':
				test_3();
				break;
			case '4':
				test_4();
				break;
			case '5':
				newTest();
				break;
			default:
				cout << endl << "\tMAIN MENU";
				cout << endl << "To start standard test enter its number(1-4).";
				cout << endl << "To create new test press 5";
				cout << endl << "To exit press \'q\'!" << endl;
				break;
		}
	}

	// cout << measurePerfomance();
	// unsigned long i = 0;
	// chrono::high_resolution_clock::time_point start, end;

	// start = end = chrono::high_resolution_clock::now();
	// for(; chrono::duration_cast<chrono::duration<double>>(end-start).count() < 1.0; 
	// 	end = chrono::high_resolution_clock::now(), ++i);

	// chrono::duration<double> elipsed =
	// 	chrono::duration_cast<chrono::duration<double>>(end-start);

	// cout << i << endl;
	// cout << elipsed.count() << " sec" << endl;


	// mt19937_64 mersenne(time(NULL));
	// unsigned long m;

	// for(int i = 0; i < 100000000; i++) {
	// 	if((m = mersenne()) < 1000000000000)
	// 		cout << mersenne() << ' ';
	// }

	// int i;
	// ifstream fin("input.txt");
	// while(!fin.eof()){
	// 	fin >> i;
	// 	cout << i << " ";
	// }
	// fin.close();

	// srand(1488);
	// vector<unsigned long> v;
	// for(int i = 0; i < 30; i++)
	// 	v.push_back(rand());
	// for(int i = 0; i < 30; i++)
	// 	cout << v[i] << " ";
	// cout << endl;

	// Tree *T = new Tree(v);

	// vector<unsigned long> nv = T->toVector();
	// for(int i = 0; i < nv.size(); i++)
	// 	cout << nv[i] << " ";
	// cout << endl;
	return 0;
}