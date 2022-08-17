/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Charts module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QDebug>
#include <QDir>
#include <QMessageBox>
#include <QValueAxis>

#include <cmath>

#include <fstream>
#include <sstream>
#include <string>
#include <numeric>
#include <map>
#include <list>
#include <tuple>
#include <vector>
#include <iterator>


using namespace std;

QT_CHARTS_USE_NAMESPACE

namespace statistics
{

template <typename U, typename T = double>
    T stddev(U begin, U end)
    {
        auto const size = distance(begin, end);

        T const sum = accumulate(begin, end, 0.0);
        T const mean = sum / size;

        vector<T> diff(size);
        transform(begin, end, diff.begin(), [mean](T x) { return x - mean; });

        T const sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        T const res = sqrt(sq_sum / size);

        return res;
    }

}


int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QChart *chart = new QChart();
    chart->legend()->hide();
    chart->createDefaultAxes();

    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    QMainWindow window;
    window.setCentralWidget(chartView);
    window.resize(1024, 768*2/3);

    string const currentDir = QDir::currentPath().toStdString() + "/../";

    map<string, double> galaxies;

    {
        ifstream file(currentDir + "gm.txt");

        for (std::string line; std::getline(file, line); )
        {
            std::istringstream in(line);

            string name;
            double fields[9];

            if (! (in >> name >> fields[0] >> fields[1] >> fields[2] >> fields[3] >> fields[4] >> fields[5] >> fields[6] >> fields[7] >> fields[8]))
                if (QMessageBox::critical(chartView, "Critical Error", "Error reading file.") == QMessageBox::Ok)
                    return -1;

            galaxies[name] = fields[8] * 1e9 * fields[0] * 2e30;
        }
    }

    map<string, tuple<list<double>, list<double>>> velocities;

    {
        ifstream file(currentDir + "grc.txt");

        for (std::string line; std::getline(file, line); )
        {
            std::istringstream in(line);

            string name;
            double fields[3];

            if (! (in >> name >> fields[0] >> fields[1] >> fields[2]))
                if (QMessageBox::critical(chartView, "Critical Error", "Error reading file.") == QMessageBox::Ok)
                    return -1;

            get<0>(velocities[name]).push_back(fields[1] * 3.08567758128e+19);
            get<1>(velocities[name]).push_back(fields[2] * 1000);
        }
    }

    static double const G = 6.671e-11;

    map<double, tuple<string, double, double, double, double, list<double>, list<double>, list<double>>> sortedgalaxies;

    for (auto i = galaxies.begin(); i != galaxies.end(); ++ i)
    {
        auto & xvelocities = get<0>(velocities[i->first]);
        auto & yvelocities = get<1>(velocities[i->first]);

        double shortestm = 0;
        double shortestw = 0;
        double shortesth = 0;
        double shorteststddev = numeric_limits<double>::max();
        list<double> shortestvelocities;

        for (double m = i->second / 15; m < i->second * 15; m += pow(10, floor(log10(i->second))) / 2)
            for (double h = 1e20 / 20; h < 1e20 * 20; h += 1e20 / 2)
            {
                list<double> ftvelocities;
                for (auto x = xvelocities.begin(); x != xvelocities.end(); ++ x)
                    ftvelocities.push_back(h / (m / (*x) + h) * sqrt(G * m / (*x)));

                double w = (* yvelocities.rbegin() - * ftvelocities.rbegin()) / * xvelocities.rbegin();

                list<double> ft2velocities;
                for (auto x = xvelocities.begin(), y = ftvelocities.begin(); x != xvelocities.end(); ++ x, ++ y)
                    ft2velocities.push_back((*y) + w * (*x));

                list<double> diff;
                for (auto y1 = yvelocities.begin(), y2 = ft2velocities.begin(); y1 != yvelocities.end(); ++ y1, ++ y2)
                    diff.push_back((*y1) - (*y2));

                double stddev = statistics::stddev(diff.begin(), diff.end());

                if (stddev < shorteststddev)
                {
                    shortestm = m;
                    shortestw = w;
                    shortesth = h;
                    shorteststddev = stddev;
                    shortestvelocities.swap(ft2velocities);
                }
            }

        sortedgalaxies[shorteststddev] = make_tuple(i->first, i->second, shortestm, shortestw, shortesth, xvelocities, yvelocities, shortestvelocities);
    }

    window.show();

    size_t count = 0;
    for (auto i = sortedgalaxies.begin(); i != sortedgalaxies.end(); ++ i, ++ count)
    {
        QLineSeries observed, theoretical;

        ostringstream out;
        out << count + 1 << "/" << sortedgalaxies.size() << ": " << get<0>(i->second) << " Tangential Velocity vs Radius (h = " << get<4>(i->second) << ", m/m_0 = " << get<2>(i->second) / get<1>(i->second) << ")";

        chart->setTitle(out.str().c_str());

        double xmax = 0, ymax = 0;

        for (auto x = get<5>(i->second).begin(), y = get<6>(i->second).begin(); x != get<5>(i->second).end(); ++ x, ++ y)
        {
            observed.append(*x, *y);

            if (*x > xmax)
                xmax = *x;
            if (*y > ymax)
                ymax = *y;
        }

        for (auto x = get<5>(i->second).begin(), y = get<7>(i->second).begin(); x != get<5>(i->second).end(); ++ x, ++ y)
        {
            theoretical.append(*x, *y);

            if (*y > ymax)
                ymax = *y;
        }

        chart->addSeries(& observed);
        chart->addSeries(& theoretical);

        QValueAxis axisX;
        axisX.setRange(0, xmax);
        axisX.setTickCount(10);
        axisX.setLabelFormat("%.2e");

        QValueAxis axisY;
        axisY.setRange(0, ymax);
        axisY.setTickCount(10);
        axisY.setLabelFormat("%.2e");

        chartView->chart()->addAxis(& axisX, Qt::AlignBottom);
        chartView->chart()->addAxis(& axisY, Qt::AlignLeft);

        if (QMessageBox::question(chartView, get<0>(i->second).c_str(), "Continue?") == QMessageBox::No)
            break;
    }

    return a.exec();
}
