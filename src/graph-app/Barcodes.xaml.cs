/// ***********************************************************
///	File	: GraphView.xaml.css
///	About	: Barcodes Visualization GUI
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Controls.Primitives;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Markup.Localizer;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using System.Xml.Linq;
using CHoleR.Helper;

namespace CHoleR
{
    /// Interaction logic for Barcodes.xaml
    public partial class Barcodes : Window
    {
        private bool barcodeMode;
        private Dictionary<int, List<KeyValuePair<Tuple<double, double>, List<int[][]>>>> PH_barcodes;
        //private Dictionary<Tuple<Point,Point>, Tuple<int, int>[]> PlotBarcodesDictionary = new Dictionary<Tuple<Point,Point>, Tuple<int, int>[]>();
        //private Point CurrentPoint;
        private double margin;
        private int mark;
        private double xmin;
        private double xmax;
        private double xstep;
        private double xframes;
        private double ymin;
        private double ymax;
        private double ystep;
        private double yframes;
        private double LastGridLineX;
        private double barcodeThickness;
        private double gridThickness;
        private int currentDimension;

        public Barcodes(bool _barcodeMode)
        {
            this.barcodeMode = _barcodeMode;
            InitializeComponent();
        }

        /// Startup actions
        private void Window_Loaded(object sender, RoutedEventArgs e)
        {            
            /// Set Window Size
            if(this.barcodeMode)
            {
                /**BarcodesWindow.Width = 640;
                BarcodesWindow.Height = 30;
                BarcodesWindow.canGraph.Width = 560;
                BarcodesWindow.canGraph.Height = 250;**/
            }
            else
            {
                BarcodesWindow.Width = 500;
                BarcodesWindow.Height = 450;
                BarcodesWindow.canGraph.Width = 400;
                BarcodesWindow.canGraph.Height = 400;
            }

            /// Set Window Position
            var desktopWorkingArea = System.Windows.SystemParameters.WorkArea;
            this.Left = desktopWorkingArea.Left;
            this.Top = desktopWorkingArea.Bottom - this.Height;
            
            /// Set Window Background
            this.Background = new SolidColorBrush( (Color)Application.Current.Resources["GraphBackground"] );
            canGraph.Background = new SolidColorBrush( (Color)Application.Current.Resources["GraphBackground"] );
            BBTitle.Background = new SolidColorBrush( (Color)Application.Current.Resources["GraphBackground"] );
            BBDimension.Background = new SolidColorBrush( (Color)Application.Current.Resources["GraphBackground"] );

            /// Set parameters
            margin = 20;
            mark = 5;
            xmin = margin + mark;
            xmax = canGraph.Width - 2*xmin;
            ymin = margin + mark + 10;
            ymax = canGraph.Height - ymin;
           
            /// Prepare dataset
            PH_barcodes = MainWindow.GraphFromLibrary.PH_barcodes;

            /// Start Drawing
            if(PH_barcodes.Count > 0)
            {
                /// Add dimension to ComboBox
                foreach (var dim in PH_barcodes)
                {
                    BBDimension.Items.Add(dim.Key);
                }

                /// Set current dimension
                int idx;
                if (PH_barcodes.Keys.Contains(1))
                {
                    currentDimension = 1;
                    idx = 1;
                }
                else
                {
                    currentDimension = PH_barcodes.Keys.Min();
                    idx = 0;
                }

                /// BARCODE MODE
                /// (*) Betti barcodes have "n" simplices less on Y-axis (n: element w/rank=0)
                if (this.barcodeMode)
                {
                    /// Draw!
                    DrawBarcodesAtDim(currentDimension);
                    BBDimension.SelectedIndex = idx;                    
                }

                /// DIAGRAM MODE
                /// (*) Persistence diagram shows births vs deaths
                else
                {
                    /// Set labels
                    BarcodesWindow.Title = "Persistence Diagram";
                    BBTitle.Text = "Persistence Diagram at dimension:";
                    
                    /// Draw!
                    DrawDiagramAtDim(currentDimension);
                    BBDimension.SelectedIndex = idx;  
                }
            }
            else
            {
                RenderText(canGraph, xmax/2 - 100, ymax/2, "No barcodes found.", Colors.Red, 13, FontWeights.UltraBold, FontStyles.Normal);
            }      
        }

        /// Draw barcodes!
        private void DrawBarcodesAtDim(int dimension)
        {
            ClearCanvas();

            var barcodes_at_dim = PH_barcodes[dimension];
            SortedSet<double> filters = new SortedSet<double>();
            int nsimplices = 0;
            foreach (var barcode in barcodes_at_dim)
            {
                var left = barcode.Key.Item1;
                var right = barcode.Key.Item2;
                if (!Double.IsInfinity(left)) filters.Add(left);
                if (!Double.IsInfinity(right)) filters.Add(right);
                nsimplices += barcode.Value.Count;
            }

            /// Set parameters
            double maxFilter = filters.Max;
            xframes = (maxFilter>=8)? maxFilter+1 : 8;
            yframes = (nsimplices>=4)? nsimplices : 4;
            xstep = (xmax - xmin)/xframes;
            ystep = (ymax - ymin)/yframes;
            barcodeThickness = (ymax - ymin) / (yframes * 2);
            gridThickness = (xmax - xmin) / (xframes * 30); ///0

            /// Render axi
            RenderAxis(xmin, xmax, xstep, xframes, ymin, ymax, ystep, yframes, mark, true, false, "t", "s", ref filters, false);

            /// Render barcodes
            RenderBarcodesA(xmin, xmax, xstep, xframes, ymin, ymax, ystep, yframes, ref barcodes_at_dim);
            //RenderBarcodesB(xmin, xmax, xstep, ymin, ymax, ystep);
            //RenderPointsCloud(xmin, xmax, xstep, ymin, ymax, ystep);

            /// Draw last gridline as red
            RenderLine(LastGridLineX, ymin, LastGridLineX, ymax, Brushes.Red, 2, false, "");
        }

        /// Axis rendering 
        private void RenderAxis(double xmin, double xmax, double xstep, double xframes, double ymin, double ymax, double ystep, double yframes, int mark, bool xlabel, bool ylabel, string xname, string yname, ref SortedSet<double> filters, bool filterAsIndex)
        {
            double xFontSize = (160/xframes <= 8)? 160/xframes : 8;
            double yFontSize = (160/yframes <= 8)? 120/xframes : 8;

            Brush lineColor = Brushes.DarkSlateGray;
            Brush gridColor = Brushes.DimGray;
            Color XYColor = Colors.DarkSlateGray;
            Color labelColor = Colors.DarkSlateGray;

            /// X axis
            RenderLine(xmin, ymax, xmax, ymax, lineColor, 1, false, "");
            RenderText(canGraph, xmax+20, ymax-10, xname, XYColor, 8, FontWeights.Bold, FontStyles.Italic);
            /// X indexes
            int f = 0;
            double err = 0.5;
            bool undertext = false;
            double maxFilterFactor = (filters.Max>0)? (filters.Max) : 0;
            double maxFilterPosition = (xstep*maxFilterFactor) + xmin;
            for (double x = xmin; x <= xmax+err; x += xstep)
            {                
                /// X mark rendering
                RenderLine(x, ymax - mark/2, x, ymax + mark/2, lineColor, 1, false, "");

                /// X label rendering
                ///TODO: if(filterAsIndex) => create new Filters[] filled with intervals
                if (xlabel)
                {
                    double xLabel = (!filterAsIndex) ? f++ : filters.ElementAt(f++);
                    FontWeight xFontWeight = (filters.Contains(xLabel))? FontWeights.DemiBold : FontWeights.ExtraLight;
                    double yleap = (undertext) ? ymax+5 :ymax;
                    undertext = !undertext;
                    RenderText(canGraph, x-5, yleap, xLabel.ToString(), labelColor, xFontSize, xFontWeight, FontStyles.Normal);
                }

                /// X grid rendering
                if ( (x >= maxFilterPosition-err) && (x <= maxFilterPosition+err) )
                {
                    /// Set last gridline
                    LastGridLineX = x;
                }
                else
                {
                    RenderLine(x, ymin, x, ymax, gridColor, gridThickness, true, "");
                }
            }
            
            /// Y axis
            RenderLine(xmin, ymin, xmin, ymax, lineColor, 1, false, "");
            RenderText(canGraph, xmin-7, ymin-20, yname, XYColor, 8, FontWeights.Bold, FontStyles.Italic);
            /// Y indexes
            double yindex = yframes;
            for (double y = ymin; y < ymax; y += ystep)
            {
                /// Y mark rendering
                RenderLine(xmin - mark/2, y, xmin + mark/2, y, lineColor, 1, false, "");
                
                /// Y label rendering
                if(ylabel)
                {
                    double yLabel = yindex--;
                    RenderText(canGraph, (xmin - mark/2) - 23, (y - yFontSize), yLabel.ToString(), labelColor, yFontSize, FontWeights.ExtraLight, FontStyles.Normal);
                }

                /// Y grid rendering
                RenderLine(xmin, y, xmax, y, gridColor, gridThickness, true, "");
            }
        }

        /// Barcodes rendering "A"
        private void RenderBarcodesA(double xmin, double xmax, double xstep, double xframes, double ymin, double ymax, double ystep, double yframes, ref List<KeyValuePair<Tuple<double, double>, List<int[][]>>> data)
        {
            /// Plot!
            int y = 0;
            foreach (var barcode in data)
            {
                foreach (var generator in barcode.Value)
                {
                    /// Determine base coordinates
                    var x1 = barcode.Key.Item1 * xstep;
                    var x2 = (!Double.IsInfinity(barcode.Key.Item2)) ? (barcode.Key.Item2 * xstep) : (xmax-xmin+10) ;
                    var y12 = ((yframes - 1) - y) * ystep;
                    y++;
                    double xpoint1 = xmin + x1;
                    double xpoint2 = xmin + x2;
                    double ypoint = ymin + y12;

                    /// Hardcode "generator" into Path.Name
                    string intervalFrom = barcode.Key.Item1.ToString();
                    string intervalTo = (barcode.Key.Item2 < Double.PositiveInfinity)? barcode.Key.Item2.ToString() : "inf";
                    string name = "_" + intervalFrom + "_" + intervalTo + "_";
                    if(currentDimension == 0)
                    {
                        name += y.ToString()+"_";
                    }
                    else
                    {
                        foreach (var element in generator)
                        {
                            foreach (var vertex in element)
                            {
                                name += vertex.ToString()+"_";
                            }                           
                        }
                    }

                    /// Render barcode for vertex "i"
                    Brush brushColor = new SolidColorBrush( MainWindow.GetSequentialColor(y, false, 0));
                    RenderLine(xpoint1, ypoint, xpoint2, ypoint, brushColor, barcodeThickness, false, name);
                    
                    /// Add line to the List of Plot Barcodes
                    /*Point point1 = new Point(xpoint1, ypoint);
                    Point point2 = new Point(xpoint2, ypoint);
                    Tuple<Point, Point> line = new Tuple<Point, Point>(point1, point2);
                    Tuple<int, int>[] check;
                    if(!PlotBarcodesDictionary.TryGetValue(line, out check))
                        PlotBarcodesDictionary.Add(line, generator);*/
                }
            }
        }

        /// Barcodes rendering "B"
        private void RenderBarcodesB(double xmin, double xmax, double xstep, double ymin, double ymax, double ystep)
        {
            #region dataset
            double[] Filters = new double[8]{1.15, 2.45, 3.37, 4.56, 5.41, 6.23, 7.13, 8.00};
            int[][] Holes = new int[8][];
            Holes[0] = new int[]{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0};
            Holes[1] = new int[]{1,1,1,1,1,1,0,1,1,1,1,1,1,1,0,1,0,0,0};
            Holes[2] = new int[]{1,0,1,1,1,1,0,1,1,1,1,1,0,1,0,1,0,0,0};
            Holes[3] = new int[]{1,0,0,1,1,1,0,1,0,1,1,1,0,1,0,0,0,0,0};
            Holes[4] = new int[]{1,0,0,1,1,0,0,1,0,1,1,0,0,1,0,0,0,0,0};
            Holes[5] = new int[]{1,0,0,1,0,0,0,1,0,1,1,0,0,1,0,0,0,0,0};
            Holes[6] = new int[]{1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0};
            Holes[7] = new int[]{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            #endregion

            int xframes = Filters.Length;
            int yframes = Holes[0].Length - 1; ///MainWindow.GraphToVisualize.VertexCount - 1;           
            xstep = (xmax - xmin)/xframes;
            ystep = (ymax - ymin)/yframes;

            Brush[] brushes = { Brushes.Red, Brushes.Blue, Brushes.DarkGoldenrod, Brushes.Green, Brushes.OrangeRed, Brushes.Purple };
            for (int i = 0; i < xframes; ++i)
            {
                for (int j = 0; j < Holes[0].Length; ++j)
                {
                    if (Holes[i][j] == 1)
                    {
                        /// Determine indexes
                        int posi = i;                   ///keep i to start from left X
                        int posj = (yframes - 1) - j; ///invert j to start from bottom Y

                        ///draw barcode for vertex "i"
                        double xpoint = xmin + (xstep*posi);
                        double ypoint = ymin + (ystep*posj);
                        GeometryGroup barcode = new GeometryGroup();
                        barcode.Children.Add(new LineGeometry(new Point(xpoint, ypoint), new Point(xpoint+xstep, ypoint)));
                        Path bar = new Path();
                        bar.StrokeThickness = 3;
                        bar.Stroke = brushes[(posj+brushes.Length) % brushes.Length];
                        bar.Data = barcode;
                        canGraph.Children.Add(bar);
                    }
                }
            }
            RenderText(canGraph, xmax/2, ymax/2, "THIS IS A TEST", Colors.DimGray, 13, FontWeights.UltraBold, FontStyles.Normal);
        }

        /// PointsCloud rendering
        private void RenderPointsCloud(double xmin, double xmax, double xstep, double ymin, double ymax, double ystep)
        {
            Brush[] brushes = { Brushes.Red, Brushes.Green, Brushes.Blue };
            Random rand = new Random();
            for (int data_set = 0; data_set < brushes.Length; data_set++)
            {
                int last_y = rand.Next((int)ymin, (int)ymax);

                PointCollection points = new PointCollection();
                for (double x = xmin; x <= xmax; x += xstep)
                {
                    last_y = rand.Next(last_y - 10, last_y + 10);
                    if (last_y < ymin) last_y = (int)ymin;
                    if (last_y > ymax) last_y = (int)ymax;
                    points.Add(new Point(x, last_y));
                }

                Polyline polyline = new Polyline();
                polyline.StrokeThickness = 1;
                polyline.Stroke = brushes[data_set];
                polyline.Points = points;

                canGraph.Children.Add(polyline);
            }
        }

        /// Draw diagram!
        private void DrawDiagramAtDim(int dimension)
        {
            ClearCanvas();

            /// Determine barcodes & max
            var barcodes_at_dim = PH_barcodes[dimension];
            /**var births = barcodes_at_dim.Select((x => x.Key.Item1)); /// .keys.births
            var deaths = barcodes_at_dim.Select((x => x.Key.Item2)).Where(d => d != double.PositiveInfinity); /// .keys.deaths (excluding INF)
            SortedSet<double> filters = new SortedSet<double>(births.Union(deaths));
            double max = filters.Max();**/
            SortedSet<double> filters = new SortedSet<double>( MainWindow.GraphFromLibrary.PH_Filtration.Keys.Select(x => (double)x));
            double maxFilter = MainWindow.GraphFromLibrary.PH_Filtration.Keys.Last();

            /// Set parameters
            xframes = (maxFilter>=4)? maxFilter+4 : 8;
            yframes = (maxFilter>=4)? maxFilter+3 : 8;
            xstep = (xmax - xmin)/xframes;
            ystep = (ymax - ymin)/yframes;
            gridThickness = (xmax - xmin) / (xframes * 90);

            /// Render axi
            RenderAxis(xmin, xmax, xstep, xframes, ymin, ymax, ystep, yframes, mark, true, true, "births", "deaths", ref filters, false);

            /// Draw bisection line
            RenderLine(xmin, ymax, xmax, ymin, Brushes.DarkRed, 1, false, "");

            /// Render barcodes
            RenderFigures(xmin, xmax, xstep, xframes, ymin, ymax, ystep, yframes, ref barcodes_at_dim);
        }

        /// Diagram figures
        private void RenderFigures(double xmin, double xmax, double xstep, double xframes, double ymin, double ymax, double ystep, double yframes, ref List<KeyValuePair<Tuple<double, double>, List<int[][]>>> data)
        {
            /// Plot!
            int y = 0;
            foreach (var barcode in data)
            {
                foreach (var generator in barcode.Value)
                {
                    y++;

                    /// Determine base coordinates
                    var birth = barcode.Key.Item1;
                    var death = barcode.Key.Item2;
                    double xpoint = xmin + birth*xstep;
                    double ypoint =  (!Double.IsInfinity(death))? (ymax - death*ystep) :  (ymin);

                    /// Hardcode "generator" into Path.Name
                    string intervalFrom = birth.ToString();
                    string intervalTo = (death < Double.PositiveInfinity)? death.ToString() : "inf";
                    string name = "_" + intervalFrom + "_" + intervalTo + "_";
                    if(currentDimension == 0)
                    {
                        name += y.ToString()+"_";
                    }
                    else
                    {
                        foreach (var element in generator)
                        {
                            foreach (var vertex in element)
                            {
                                name += vertex.ToString()+"_";
                            }                           
                        }
                    }

                    /// Render barcode for vertex "i"
                    Brush brushColor = new SolidColorBrush( MainWindow.GetSequentialColor(y, false, 0));
                    double size = 6;///xstep/6;
                    if(currentDimension != 0)
                    {
                        PointCollection myPointCollection = new PointCollection();
                        myPointCollection.Add(new Point(xpoint - size, ypoint + size));
                        myPointCollection.Add(new Point(xpoint + 2*size - size, ypoint + size));
                        myPointCollection.Add(new Point(xpoint, ypoint - 2*size + size));
                        RenderPolygon(myPointCollection, brushColor, 0, name);
                    }
                    else
                    {
                        RenderEllipse(new Point(xpoint - size/2, ypoint - size/2), brushColor, size, 0, name);
                    }
                }
            }
        }

        #region TOOLS
        /// Render line
        private void RenderLine(double xpoint1, double ypoint1, double xpoint2, double ypoint2, Brush colour, double thickness, bool dashed, string name)
        {
            /// Draw line
            GeometryGroup figure = new GeometryGroup();
            Point point1 = new Point(xpoint1, ypoint1);
            Point point2 = new Point(xpoint2, ypoint2);
            figure.Children.Add(new LineGeometry(point1, point2));
            Path line = new Path();
            line.StrokeThickness = thickness;
            if(dashed) line.StrokeDashArray = new DoubleCollection() { 2, 4 };  ///{length, space}
            line.Stroke = colour;
            line.Data = figure;
            line.Name = name.ToString();
            canGraph.Children.Add(line);
        }

        /// Render text
        private void RenderText(Canvas canvas, double x, double y, string text, Color color, double size, FontWeight weight, FontStyle style)
        {
            Point point = new Point(x, y);
            
            var label = new Label
            {
                Content = text,
                Margin = new Thickness(point.X, point.Y, 0, 0),
                HorizontalContentAlignment = HorizontalAlignment.Center,
                VerticalContentAlignment = VerticalAlignment.Center,
                Background = new SolidColorBrush(Colors.Transparent),
                Foreground = new SolidColorBrush(color),
                FontSize = size,
                FontWeight = weight,
                FontStyle = style
            };
            canvas.Children.Add(label);
        }

        /// Render polygon
        private void RenderPolygon(PointCollection myPointCollection, Brush brushColor, double thickness, string name)
        {
            Polygon myPolygon = new Polygon();
            //myPolygon.Name = "simplex_" + myPointCollection.GetHashCode().ToString();
            myPolygon.Name = name;
            myPolygon.Points = myPointCollection;
            myPolygon.Opacity = 0.8;
            myPolygon.Fill = brushColor;
            myPolygon.Stroke = Brushes.Transparent;
            myPolygon.StrokeThickness = thickness;
                    
            canGraph.Children.Add(myPolygon);
        }

        ///  Render ellipse
        private void RenderEllipse(Point myPoint, Brush brushColor, double diameter, double thickness, string name)
        {
            Ellipse ellipse = new Ellipse();
            ellipse.Name = name;
            ellipse.Fill = brushColor;
            ellipse.Stroke = Brushes.Transparent;
            ellipse.StrokeThickness = thickness;
            Canvas.SetZIndex(ellipse, (int)thickness);
            ellipse.Height = diameter;
            ellipse.Width = diameter;
            ellipse.Margin = new Thickness(myPoint.X, myPoint.Y, 0, 0);

            canGraph.Children.Add(ellipse);
        }

        /// Clear canvas (selectively)
        private void ClearCanvas()
        {
            /// Clear canGraph.Children
            for (int index = canGraph.Children.Count - 1; index >= 0; index--)
            {
                if ( !(canGraph.Children[index] is TextBox) && !(canGraph.Children[index] is ComboBox))
                {
                    ///canGraph.Children.Remove((UIElement)children);
                    canGraph.Children.RemoveAt(index);
                }
            }
        }

        #endregion

        #region EVENTS

        /// Zoom canvas
        private void sliZoom_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            /// Make sure the control's are all ready
            if (!IsInitialized) return;

            /// Display the zoom factor as a percentage
            lblZoom.Content = sliZoom.Value + "%";

            /// Get the scale factor as a fraction 0.25 - 2.00
            double scale = (double)(sliZoom.Value / 100.0);

            /// Scale the graph.
            canGraph.LayoutTransform = new ScaleTransform(scale, scale);
        }

        ///  Zoom on window resize
        private void BarcodesWindow_SizeChanged(object sender, SizeChangedEventArgs e)
        {
            switch (this.WindowState) 
            { 
                case WindowState.Maximized:
                    sliZoom.Value = 200;
                    break; 
                case WindowState.Minimized: 
                    sliZoom.Value = 100;
                    break; 
                case WindowState.Normal: 
                    sliZoom.Value = 100;
                    break; 
            } 
        }

        /// On mouse click
        private ToolTip CanvasToolTip = new ToolTip(){IsOpen = false};
        private void Canvas_Click(object sender, MouseEventArgs e)
        {
            /// Based on clicked path
            if ( (e.OriginalSource is Path) || (e.OriginalSource is Polygon) || (e.OriginalSource is Ellipse))
            {
                if (e.OriginalSource is Path)
                {
                    Path simplex = (Path) e.OriginalSource;

                    if(simplex.Name.StartsWith("_"))
                    {
                        ///Get data from path
                        var pathStroke = simplex.Stroke;
                        Color color = (pathStroke as SolidColorBrush).Color;
                        var pathX = simplex.Data.Bounds.X;
                        var pathXWidth = simplex.Data.Bounds.Width;
                        var pathY = simplex.Data.Bounds.Y;
                        var pathYHeight = simplex.Data.Bounds.Height;
                        var pathName = simplex.Name;
                        string[] generatorS = pathName.Split('_').ToArray();
                        string intervalFrom = generatorS[1];
                        string intervalTo = generatorS[2];
                        var generators = generatorS.Skip(3).Take(generatorS.Count()-4).ToList();
                
                        /// Interval
                        string closure = (intervalTo != "inf") ? "]" : ")";
                        string text = "[" + intervalFrom + "," + intervalTo + closure;

                        ///Show ToolTip
                        CanvasToolTip.IsOpen = false;
                        canGraph.ToolTip = CanvasToolTip;
                        CanvasToolTip.Content = text;
                        CanvasToolTip.Background = pathStroke;
                        CanvasToolTip.Foreground = Brushes.WhiteSmoke;
                        CanvasToolTip.IsOpen = true;
                        ToolTipService.SetShowDuration(canGraph, 2000);
                
                        /// Highligh simplices
                        MainWindow.RemoveCustomHighlights();
                        int k = currentDimension + 1;
                        for (int i = 0; i < generators.Count; i += k)
                        {
                            int[] elements = generators.Skip(i).Take(k).Select(int.Parse).ToArray();
                            MainWindow.HighlightPath(ref elements, 2, 1, 1, true, false, color);
                        }
                    }
                }
                else if ( e.OriginalSource is Polygon)
                {
                    Polygon simplex = (Polygon) e.OriginalSource;

                    if(simplex.Name.StartsWith("_"))
                    {
                        ///Get data from polygon
                        var polygonFill = simplex.Fill;
                        Color color = (polygonFill as SolidColorBrush).Color;

                        var pathName = simplex.Name;
                        string[] generatorS = pathName.Split('_').ToArray();
                        string intervalFrom = generatorS[1];
                        string intervalTo = generatorS[2];
                        var generators = generatorS.Skip(3).Take(generatorS.Count()-4).ToList();
                
                        /// Interval
                        string closure = (intervalTo != "inf") ? "]" : ")";
                        string text = "[" + intervalFrom + "," + intervalTo + closure;

                        ///Show ToolTip
                        CanvasToolTip.IsOpen = false;
                        canGraph.ToolTip = CanvasToolTip;
                        CanvasToolTip.Content = text;
                        CanvasToolTip.Background = polygonFill;
                        CanvasToolTip.Foreground = Brushes.WhiteSmoke;
                        CanvasToolTip.IsOpen = true;
                        ToolTipService.SetShowDuration(canGraph, 2000);
                
                        /// Highligh simplices
                        MainWindow.RemoveCustomHighlights();
                        int k = currentDimension + 1;
                        for (int i = 0; i < generators.Count; i += k)
                        {
                            int[] elements = generators.Skip(i).Take(k).Select(int.Parse).ToArray();
                            MainWindow.HighlightPath(ref elements, 2, 1, 1, true, false, color);
                        }
                    }
                }
                else
                {
                    Ellipse simplex = (Ellipse) e.OriginalSource;

                    if(simplex.Name.StartsWith("_"))
                    {
                        ///Get data from polygon
                        var polygonFill = simplex.Fill;
                        Color color = (polygonFill as SolidColorBrush).Color;

                        var pathName = simplex.Name;
                        string[] generatorS = pathName.Split('_').ToArray();
                        string intervalFrom = generatorS[1];
                        string intervalTo = generatorS[2];
                        var generators = generatorS.Skip(3).Take(generatorS.Count()-4).ToList();
                
                        /// Interval
                        string closure = (intervalTo != "inf") ? "]" : ")";
                        string text = "[" + intervalFrom + "," + intervalTo + closure;

                        ///Show ToolTip
                        CanvasToolTip.IsOpen = false;
                        canGraph.ToolTip = CanvasToolTip;
                        CanvasToolTip.Content = text;
                        CanvasToolTip.Background = polygonFill;
                        CanvasToolTip.Foreground = Brushes.WhiteSmoke;
                        CanvasToolTip.IsOpen = true;
                        ToolTipService.SetShowDuration(canGraph, 2000);
                
                        /// Highligh simplices
                        MainWindow.RemoveCustomHighlights();
                        int k = currentDimension + 1;
                        for (int i = 0; i < generators.Count; i += k)
                        {
                            int[] elements = generators.Skip(i).Take(k).Select(int.Parse).ToArray();
                            MainWindow.HighlightPath(ref elements, 2, 1, 1, true, false, color);
                        }
                    }
                }
            }

            /// Based on clicked point
            /*Point _currentPoint = e.GetPosition(this);
            CurrentPoint = _currentPoint;
            double pointX = _currentPoint.X;
            double pointY = _currentPoint.Y;
            Console.WriteLine("Left click at (" + pointX + "," + pointY + ")");

            foreach (var barcode in PlotBarcodesDictionary)
            {
                var line = barcode.Key;
                double yErr = barcodeThickness / 2;
                double offset = 1.51; ///TODO
                double x1 = line.Item1.X;
                double x2 = line.Item2.X;
                double y12 = line.Item1.Y;

                if ((x1 < xmax) && (x2 < xmax + xmin) && (y12 < ymax)
                    && (pointX >= x1 && pointX <= x2)
                    && (pointY >= y12 - yErr && pointY <= y12 + yErr))
                {
                    double start = x1 / xstep;
                    double end = (x2 < (xmax - xmin + 10)) ? x2 / xstep : Double.PositiveInfinity;
                    ///string text = "[" + Math.Truncate(start*100)/100 + " -> " + Math.Truncate(end*100)/100 + "]";                 
                    string text = barcode.Value.Length + "-simplex";

                    ///Show TT
                    canGraph.ToolTip = CanvasToolTip;
                    CanvasToolTip.Content = text;
                    CanvasToolTip.IsOpen = true;
                    ToolTipService.SetShowDuration(canGraph, 2000);

                    ///Hide TT
                    //canGraph.ToolTip = null;
                    CanvasToolTip.IsOpen = false;
                    break;
                }
            }*/
        }

        /// On selection changed
        private void BBDimension_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (BBDimension.SelectedItem != null)
            {
                currentDimension = int.Parse(BBDimension.SelectedItem.ToString());
            }

            if(this.barcodeMode)
                DrawBarcodesAtDim(currentDimension);
            else
                DrawDiagramAtDim(currentDimension);
        }

        /// On windows exit
        private void BarcodesWindow_Closing(object sender, CancelEventArgs e)
        {
            ///BarcodesWindow.Dispatcher.Thread.Interrupt();
            CanvasToolTip.IsOpen = false;
        }

        /// Entropy button
        private void BtnEntropy_Click(object sender, RoutedEventArgs e)
        {
            /// Compute persistent entropy: 
            /// "sum all barcodes (if inf replace with max+1) over all dimensions, then divide by number of dimensions??"
            /*int ndimensions = PH_barcodes.Count;
            int nbarcodes = 0;
            foreach (var barcode in PH_barcodes.Values)
            {
                nbarcodes += barcode.Count;
            }*/

            /// Determine maximum end among all barcodes        
            double bMax = 0;
            foreach (var barcode in PH_barcodes.Values)
            {
                /*var startPoints = barcode.Keys;
                double iMaxStart =  startPoints.Max(y => y.Item1);

                var endPoints = barcode.Keys.Where(y => y.Item2 != double.PositiveInfinity); ///exclude INF
                double iMaxEnd = (endPoints.Any()) ? endPoints.Max(y => y.Item2) : 0;

                ///double iMax = (intervals.Count() > 0)? dimension.Keys.Max(y => y.Item2) : 0; /// Maximum among EndPoints
                double iMax = (iMaxStart > iMaxEnd)? iMaxStart : iMaxEnd;
                bMax = (iMax > bMax)? iMax : bMax;*/
            
                var startPoints = barcode.Select(x => x.Key).ToArray(); /// .keys
                double iMaxStart =  startPoints.Max(y => y.Item1);

                var endPoints = barcode.Select(x => x.Key).Where(y => y.Item2 != double.PositiveInfinity); /// .keys excluding INF
                double iMaxEnd = (endPoints.Any()) ? endPoints.Max(y => y.Item2) : 0;

                ///double iMax = (intervals.Count() > 0)? dimension.Keys.Max(y => y.Item2) : 0; /// Maximum among EndPoints
                double iMax = (iMaxStart > iMaxEnd)? iMaxStart : iMaxEnd;
                bMax = (iMax > bMax)? iMax : bMax;
                
            }

            /// Determine sum of all lengths
            double bSum = 0;
            double bLong = 0;
            foreach (var barcode in PH_barcodes.Values)
            {
                foreach (var interval in barcode)
                {
                    double bStart = interval.Key.Item1;
                    double bEnd = (interval.Key.Item2 == double.PositiveInfinity)? bMax+1 : interval.Key.Item2;
                    double bLength = bEnd - bStart;
                    bSum += bLength;
                    bLong = (bLength > bLong) ? bLength : bLong;
                }
            }

            /// Determine entropies
            Dictionary<int, double> entropies = PH_barcodes.Keys.ToArray().ToDictionary(v => v, v => 0.00);
            string result = "";
            foreach (var barcodesAtDim in PH_barcodes)
            {
                int dim = barcodesAtDim.Key;
                foreach (var barcode in barcodesAtDim.Value)
                {
                    double bstart = barcode.Key.Item1;
                    double bend = (barcode.Key.Item2 == double.PositiveInfinity)? bMax+1 : barcode.Key.Item2;
                    double blength = bend - bstart;
                    double p = ( blength / bSum ) * ( Math.Log(blength/bSum) );
                    entropies[dim] += -p;
                }
                result +=  "P.E. (" + dim + ") = " +  entropies[dim].ToString() +"\n";
            }
            double PE = entropies.Sum(x => x.Value);
            result += "\n Standard P.E. = " + PE + "\n";

            double NPE = PE / Math.Log(bLong);
            result += "\n Normalized P.E. = " + NPE + "\n";

            /// Display results (more uniform then less entropy)
            MessageBoxResult messageBoxResult = System.Windows.MessageBox.Show(
                "\n" + result + "\n", "Persistent entropy",
                System.Windows.MessageBoxButton.OK, MessageBoxImage.Information);

            /// Graph 20 0.20
            /// P.E at Homology: 0= 1.936366655135545
            /// P.E at Homology: 1= 1.359085121923069
            /// P.E.=  3.2954517770586143
            
            /// Graph stella 8
            /// P.E at Homology: 0= 1.7498125725957974
            /// P.E at Homology: 1= 0.1037164860181224
            /// P.E.=  1.8535290586139197
        }
        
        #endregion
    }
}
