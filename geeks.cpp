// C++ program for Implementation of
// Reverse Cuthill Mckee Algorithm

#include <bits/stdc++.h>
using namespace std;

vector<double> globalDegree;

int findIndex(vector<pair<int, double> > a, int x)
{
    for (int i = 0; i < a.size(); i++)
        if (a[i].first == x)
            return i;
    return -1;
}

bool compareDegree(int i, int j)
{
    return ::globalDegree[i] < ::globalDegree[j];
}

template <typename T>
ostream& operator<<(ostream& out, vector<T> const& v)
{
    for (int i = 0; i < v.size(); i++)
        out << v[i] << ' ';
    return out;
}

class ReorderingSSM {
private:
    //vector<vector<double> > _matrix;
    int nodesize;
    int *matrixRowptr, *matrixColind;
    double *matrixValues;

public:
    // Constructor and Destructor
    ReorderingSSM(int x, double *a, int*z, int*t) : nodesize(x), matrixValues(a), matrixRowptr(z), matrixColind(t)
    {}

    ~ReorderingSSM() {}

    // class methods

    // Function to generate degree of all the nodes
    vector<double> degreeGenerator()
    {

        vector<double> degrees;
        int elmCount;

        for (int i = 0; i < nodesize; i++) {
            double count = 0;
            elmCount = matrixRowptr[i+1] - matrixRowptr[i];
            double *pos = matrixValues + matrixRowptr[i];
            for (int j = 0; j < elmCount; pos++, j++) {
                count += *pos;
            }

            degrees.push_back(count);
        }

        return degrees;
    }

    // Implementation of Cuthill-Mckee algorithm
    vector<int> CuthillMckee()
    {
        vector<double> degrees = degreeGenerator();
        cout << "degreeGenerator finished." << endl << flush;

        ::globalDegree = degrees;

        queue<int> Q;
        vector<int> R;
        vector<pair<int, double> > notVisited;

        for (int i = 0; i < degrees.size(); i++)
            notVisited.push_back(make_pair(i, degrees[i]));

        // Vector notVisited helps in running BFS
        // even when there are dijoind graphs
        while (notVisited.size()) {

            int minNodeIndex = 0;

            for (int i = 0; i < notVisited.size(); i++)
                if (notVisited[i].second < notVisited[minNodeIndex].second)
                    minNodeIndex = i;

            Q.push(notVisited[minNodeIndex].first);

            notVisited.erase(notVisited.begin()
                             + findIndex(notVisited,
                                         notVisited[Q.front()].first));

            // Simple BFS
            cout << "before bfs." << endl << flush;
            while (!Q.empty()) {

                vector<int> toSort;

                int front = Q.front();
                int neighbourCount = matrixRowptr[front+1] - matrixRowptr[front];
                double *pos = matrixValues + matrixRowptr[front];
                int *colIndPos = matrixColind + matrixRowptr[front];

                for (int i = 0; i < neighbourCount; colIndPos++, pos++, i++) {
                    if (*colIndPos != front && findIndex(notVisited, *colIndPos) != -1) {
                        toSort.push_back(*colIndPos);
                        notVisited.erase(notVisited.begin()
                                         + findIndex(notVisited, *colIndPos));
                    }
                }
                cout << "after bfs." << endl << flush;

                sort(toSort.begin(), toSort.end(), compareDegree);

                for (int i = 0; i < toSort.size(); i++)
                    Q.push(toSort[i]);

                cout << "after sort." << endl << flush;

                R.push_back(Q.front());
                Q.pop();
            }
        }

        return R;
    }

    // Implementation of reverse Cuthill-Mckee algorithm
    vector<int> ReverseCuthillMckee()
    {

        vector<int> cuthill = CuthillMckee();
        cout << "algo finished." << endl << flush;

        int n = cuthill.size();
        cout << "cuthill size: " << n << endl << flush;

        if (n % 2 == 0)
            n -= 1;

        n = n / 2;

        for (int i = 0; i <= n; i++) {
            int j = cuthill[cuthill.size() - 1 - i];
            cuthill[cuthill.size() - 1 - i] = cuthill[i];
            cuthill[i] = j;
        }

        return cuthill;
    }
};
