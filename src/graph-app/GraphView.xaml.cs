/// ***********************************************************
///	File	: GraphView.xaml.css
///	About	: Graph Visualization GUI
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using System.Windows.Threading;
using CHoleR.Helper;
using GRAPH;
using QuickGraph;
using GraphSharp;
using GraphSharp.Algorithms;
using GraphSharp.Controls;

namespace CHoleR
{
        public partial class GraphView : Window
        {
            #region CLASS VARIABLES

            public GraphManager GraphFromLibrary;

            private BidirectionalWeightedGraph<object, WeightedEdge<object>> _graphToVisualize;
            public BidirectionalWeightedGraph<object, WeightedEdge<object>> GraphToVisualize
            {
                get { return this._graphToVisualize; }
                set
                {
                    if (!Equals(value, this._graphToVisualize))
                    {
                        this._graphToVisualize = value;
                        this.RaisePropChanged("GraphToVisualize");
                    }
                }
            }

            public event PropertyChangedEventHandler PropertyChanged;
            public void RaisePropChanged(string name)
            {
                var eh = this.PropertyChanged;
                if (eh != null)
                {
                    eh(this, new PropertyChangedEventArgs(name));
                }
            }

            public volatile bool MouseDownOnGrid;
            
            public List<GraphCycle> GraphCycles {get;set;}
            public class GraphCycle
            {
                public double Weight {get;set;}
                public int Edges {get;set;}
                public string Cycle {get;set;}
            };

            public List<GraphClique> GraphCliques {get;set;}
            public class GraphClique
            {
                public int Number {get;set;}
                public string Vertices {get;set;}
            };

            public int ShowingNow = 0;
            
            #endregion
            
            #region INITIALIZATION
            public GraphView(ref GraphManager G, ref BidirectionalWeightedGraph<object, WeightedEdge<object>> graph)
            {
                /// ASSIGN graph to viewer
                GraphFromLibrary = G;
                GraphToVisualize = graph;

                InitializeComponent();

                /// Customize layout
                SetCustomLayout(GraphFromLibrary.V, GraphFromLibrary.E);
                MainWindow.graphLayout = graphLayout;

                /// Set background
                //this.Background = new SolidColorBrush( (Color)Application.Current.Resources["GraphBackground"] );
                //GraphZoom.Background = new SolidColorBrush( (Color)Application.Current.Resources["GraphBackground"] );
                ImageBrush Brushes_Background = new ImageBrush();
                Brushes_Background.ImageSource = new BitmapImage( new Uri("pack://application:,,,/"+ Assembly.GetExecutingAssembly().GetName().Name + ";component/Resources/graph_bg.jpg") );
                GraphZoom.Background = Brushes_Background;
                //graphLayout.Background = new SolidColorBrush( (Color)Application.Current.Resources["GraphBackground"] );
                graphLayout.Background = Brushes.Transparent;
                
                /// Clear HighlightedEdges
                MainWindow.HighlightedEdges.Clear();
                MainWindow.HighlightedVertices.Clear();

                /// SET graph information
                if(GraphFromLibrary.source != "")
                {
                    string gsource = GraphFromLibrary.source.Remove(0, 6); // remove "input_"
                    LabelInfo_1B.Content = gsource.Remove(gsource.Length - 4); // remove ".txt

                    LabelInfo_2B.Content = GraphFromLibrary.V;
                    LabelInfo_3B.Content = GraphFromLibrary.E;
                    LabelInfo_4B.Content = Math.Round(GraphFromLibrary.Density, 3);
                    LabelInfo_5B.Content = GraphFromLibrary.Degeneracy;
                    LabelInfo_6B.Content = GraphFromLibrary.Cyclomatic;
                    LabelInfo_7B.Content = GraphFromLibrary.Components;
                }

                /// SHOW rendering algorithms
                ///{ "CompoundFDP", "BoundedFR", "FR", "LinLog", "KK", "Circular", "Tree", "ISOM", "EfficientSugiyama" };
                IEnumerable<string> RenderingAlgorithms = new string[] { "CompoundFDP", "BoundedFR", "LinLog", "KK", "Circular", "Tree", "ISOM" };
                foreach (string algo in RenderingAlgorithms)
                {
                    this.ComboAlgorithm.Items.Add(algo);
                }
                this.ComboAlgorithm.SelectedIndex = 3;

                /// SHOW cycles and cliques (if available)
                this.BtnCycles.Visibility = (GraphFromLibrary.minimumCycleBasis != null && GraphFromLibrary.minimumCycleBasis.Length > 0)? System.Windows.Visibility.Visible : System.Windows.Visibility.Hidden;
                this.BtnCliques.Visibility = (GraphFromLibrary.maximalCliqueSet != null && GraphFromLibrary.maximalCliqueSet.Length > 0)? System.Windows.Visibility.Visible : System.Windows.Visibility.Hidden;
                
                /// SHOW persistent homology (ifavailable)
                this.BtnBetti.Visibility = (GraphFromLibrary.PH_barcodes != null && GraphFromLibrary.PH_barcodes.Count > 0)? System.Windows.Visibility.Visible : System.Windows.Visibility.Hidden;
                this.BtnDiagram.Visibility = (GraphFromLibrary.PH_barcodes != null && GraphFromLibrary.PH_barcodes.Count > 0)? System.Windows.Visibility.Visible : System.Windows.Visibility.Hidden;
                this.BtnBottleneck.Visibility = (GraphFromLibrary.PH_barcodes != null && GraphFromLibrary.PH_barcodes.Count > 0)? System.Windows.Visibility.Visible : System.Windows.Visibility.Hidden;
                this.BtnFiltration.Visibility = (GraphFromLibrary.PH_barcodes != null && GraphFromLibrary.PH_barcodes.Count > 0)? System.Windows.Visibility.Visible : System.Windows.Visibility.Hidden;

                /// Windows size depends on displayed vertices
                //this.WindowState = (G.V<10)? WindowState.Normal : WindowState.Maximized;
                this.WindowState = WindowState.Maximized;
            }
        
            #endregion

            #region BUTTON ACTIONS

            /// Change view
            private void BtnChangeView_Click(object sender, RoutedEventArgs e)
            { 
                /// Update layout
                UpdateGraphLayout();

                /// Remove filtrations
                MainWindow.RemoveSimplices();

                /// Reset GUI
                SliderFilters.Value = SliderFilters.Minimum;
                FiltrationLabel.Content = "";
            }

            /// Button show cycles
            private void BtnShowCycles_Click(object sender, RoutedEventArgs e)
            {
                if(ShowingNow == 1)
                {
                    ShowingNow = 0;
                    ///BtnCycles.Background  = new SolidColorBrush(Colors.DimGray);
                    BtnCycles.IsEnabled = true;
                    ///BtnCycles.Content = "Show cycles";
                    DataGrid1.Visibility = System.Windows.Visibility.Hidden; 
                    BtnShowAll.Visibility = System.Windows.Visibility.Hidden; 
                }
                else
                {
                    ShowingNow = 1;
                    ///BtnCycles.Background  = new SolidColorBrush(Colors.Green);
                    BtnCycles.IsEnabled = false;
                    ShowCyclesList();
                    ///BtnCycles.Content = "Hide cycles";
                    DataGrid1.Visibility = System.Windows.Visibility.Visible; 
                    BtnShowAll.Visibility = System.Windows.Visibility.Visible;

                    ///BtnCliques.Background  = new SolidColorBrush(Colors.DimGray);
                    BtnCliques.IsEnabled = true;
                    ///BtnCliques.Content = "Show cliques";
                }
            }

            /// Button show cliques
            private void BtnShowCliques_Click(object sender, RoutedEventArgs e)
            {
                 if(ShowingNow == 2)
                {
                    ShowingNow = 0;
                    ///BtnCliques.Background  = new SolidColorBrush(Colors.DimGray);
                    BtnCliques.IsEnabled = true;
                    ///BtnCliques.Content = "Show cliques";
                    DataGrid1.Visibility = System.Windows.Visibility.Hidden; 
                    BtnShowAll.Visibility = System.Windows.Visibility.Hidden; 
                }
                else
                {
                    ShowingNow = 2;
                    ///BtnCliques.Background  = new SolidColorBrush(Colors.Green);
                    BtnCliques.IsEnabled = false;
                    ShowCliquesList();
                    ///BtnCliques.Content = "Hide cliques";
                    DataGrid1.Visibility = System.Windows.Visibility.Visible; 
                    BtnShowAll.Visibility = System.Windows.Visibility.Visible;

                    ///BtnCycles.Background  = new SolidColorBrush(Colors.DimGray);
                    BtnCycles.IsEnabled = true;
                    ///BtnCycles.Content = "Show cycles";
                }
            }

            /// Show all cycles/cliques
            private void BtnShowAll_Click(object sender, RoutedEventArgs e)
            {
                DataGrid1.IsEnabled = false;
                BtnShowAll.IsEnabled = false;

                if(ShowingNow == 1)
                {
                    foreach (var item in GraphFromLibrary.minimumCycleBasis.Select((value, idx) => new { value, idx }))
                    {
                        var cycle = item.value;
                        var index = item.idx;
                        MainWindow.HighlightPath(ref cycle, 0, GraphFromLibrary.minimumCycleBasis.Length, 0, true, true, Colors.Transparent);
                    }
                }

                if (ShowingNow == 2)
                {
                    BtnCliques.IsEnabled = false;
                    foreach (var item in GraphFromLibrary.maximalCliqueSet.Select((value, idx) => new { value, idx }))
                    {
                        var clique = item.value;
                        var index = item.idx;
                        MainWindow.HighlightPath(ref clique, 2, GraphFromLibrary.maximalCliqueSet.Length, 0, true, true, Colors.Transparent);
                    }
                    BtnCliques.IsEnabled = true;
                }

                BtnShowAll.IsEnabled = true;
                DataGrid1.IsEnabled = true;
            }

            /// Show Betti numbers
            private void BtnBettiNumbers_Click(object sender, RoutedEventArgs e)
            {
                ///Disable eventually active Cycles & Cliques
                ShowingNow = 0;
                ///BtnCycles.Background  = new SolidColorBrush(Colors.DimGray);
                BtnCycles.IsEnabled = true;
                ///BtnCliques.Background  = new SolidColorBrush(Colors.DimGray);
                BtnCliques.IsEnabled = true;
                ///BtnCycles.Content = "Show cycles";
                ///BtnCliques.Content = "Show cliques";
                DataGrid1.Visibility = System.Windows.Visibility.Hidden; 
                BtnShowAll.Visibility = System.Windows.Visibility.Hidden; 
                
                ///Show betti barcodes in a new thread
                var BarcodesWindow = new Barcodes(true);
                try
                {
                    /// SHOW() WITH "SEND" (...<Background<Input<Loaded<Render<DataBind<Normal<Send) PRIORITY              
                    Task.Run(()=> { Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Send,
                        new Action(() =>
                        {
                            BarcodesWindow.Show();
                            ///BarcodesWindow.ShowDialog();
                        }));
                    });
                }
                catch (ThreadInterruptedException ex)
                {
                    Console.WriteLine("Thread interrupted! (Barcodes) - " + ex);
                }
            }

            /// Show Persistence diagram
            private void BtnDiagram_Click(object sender, RoutedEventArgs e)
            {
                ///Disable eventually active Cycles & Cliques
                ShowingNow = 0;
                ///BtnCycles.Background  = new SolidColorBrush(Colors.DimGray);
                BtnCycles.IsEnabled = true;
                ///BtnCliques.Background  = new SolidColorBrush(Colors.DimGray);
                BtnCliques.IsEnabled = true;
                ///BtnCycles.Content = "Show cycles";
                ///BtnCliques.Content = "Show cliques";
                DataGrid1.Visibility = System.Windows.Visibility.Hidden; 
                BtnShowAll.Visibility = System.Windows.Visibility.Hidden; 
                
                ///Show betti barcodes in a new thread
                var BarcodesWindow = new Barcodes(false);
                try
                {
                    /// SHOW() WITH "SEND" (...<Background<Input<Loaded<Render<DataBind<Normal<Send) PRIORITY              
                    Task.Run(()=> { Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Send,
                        new Action(() =>
                        {
                            ///BarcodesWindow.Show();
                            BarcodesWindow.ShowDialog();
                        }));
                    });
                }
                catch (ThreadInterruptedException ex)
                {
                    Console.WriteLine("Thread interrupted! (Persistence Diagram) - " + ex);
                }
            }

            /// Display animated filtration
            private void BtnFilters_Click(object sender, RoutedEventArgs e)
            {
                ///Disable eventually active Cycles & Cliques
                ShowingNow = 0;
                ///BtnCycles.Background  = new SolidColorBrush(Colors.DimGray);
                BtnCycles.IsEnabled = true;
                ///BtnCliques.Background  = new SolidColorBrush(Colors.DimGray);
                BtnCliques.IsEnabled = true;
                ///BtnCycles.Content = "Show cycles";
                ///BtnCliques.Content = "Show cliques";
                DataGrid1.Visibility = System.Windows.Visibility.Hidden; 
                BtnShowAll.Visibility = System.Windows.Visibility.Hidden;
                BorderFilters.Margin = new Thickness{ Left = (this.Width/2), Right = 0, Top = 0, Bottom = 0};

                /// Set slider
                int maxFilters = GraphFromLibrary.PH_Filtration.Count;
                SliderFilters.Maximum = maxFilters;
                SliderFilters.Value = SliderFilters.Minimum;

                /// Show slider
                BorderFilters.Visibility = Visibility.Visible;

                /// Trigger filtration
                ///FiltrationLabel.Content = 
                MainWindow.AnimatedFiltration(0, 0);
            }

            /// Display bottleneck diagram
            private void BtnBottleneck_Click(object sender, RoutedEventArgs e)
            {
                ///Show bottleneck distance in a new thread
                var BottleneckWindow = new Bottleneck();
                try
                {
                    /// SHOW() WITH "SEND" (...<Background<Input<Loaded<Render<DataBind<Normal<Send) PRIORITY              
                    Task.Run(()=> { Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Send,
                        new Action(() =>
                        {
                            ///BottleneckWindow.Show();
                            BottleneckWindow.ShowDialog();
                        }));
                    });
                }
                catch (ThreadInterruptedException ex)
                {
                    Console.WriteLine("Thread interrupted! (Bottleneck) - " + ex);
                }
            }

            /// Reset graph visualization
            private void BtnReset_Click(object sender, RoutedEventArgs e)
            {
                /// Restore hidden vertices&edges
                if(graphLayout.Graph.HiddenVertices.Any())
                {
                    foreach (var vertex in graphLayout.Graph.HiddenVertices.ToList())
                    {
                        GraphToVisualize.UnhideVertexAndEdges(vertex);
                        MainWindow.HighlightedVertices.Clear();
                        MainWindow.HighlightedEdges.Clear();
                    }
                }

                /// Restore colors of vertices
                
                foreach (TaggedVertex vertex in GraphToVisualize.Vertices)
                {
                    vertex.Color = (Color) Application.Current.Resources["VertexNormal"];

                    Color vertexEnabledColor = MainWindow.HexToColor("#FFFFFF");
                    VertexControl VC = graphLayout.GetVertexControl(vertex);
                    VC.Foreground =  new SolidColorBrush(vertexEnabledColor);
                    VC.BorderBrush =  new SolidColorBrush(vertexEnabledColor);
                }
                
                /// Restore colors of edges
                foreach (TaggedEdge edge in GraphToVisualize.Edges)
                {
                    edge.Color = (Color) Application.Current.Resources["EdgeNormal"];
                }

                /// Remove highlighted filtrations
                MainWindow.RemoveSimplices();

                ///Reset GUI
                ResetGUI();
            }

            /// Reset GUI elements
            private void ResetGUI()
            {
                /// Restore buttons
                DataGrid1.Visibility = Visibility.Hidden;
                BtnShowAll.Visibility = Visibility.Hidden; 
                
                
                BtnCycles.IsEnabled = true;
                //BtnCycles.Content = "Minimum cycles";
                //BtnCycles.Background  = new SolidColorBrush(Colors.DimGray);
                
                BtnCliques.IsEnabled = true;
                //BtnCliques.Content = "Maximal cliques";
                //BtnCliques.Background  = new SolidColorBrush(Colors.DimGray);

                BtnBetti.IsEnabled = true;
                BtnFiltration.IsEnabled = true;
                BtnBottleneck.IsEnabled = true;
                BorderFilters.Visibility = Visibility.Hidden;

                ShowingNow = 0;
            }

            #endregion

            #region MAIN FUNCTIONS

            private void UpdateGraphLayout()
            {
                ///{ "Automatic", "Compound", "Simple" };
                graphLayout.LayoutMode = LayoutMode.Automatic;

                ///{ "BoundedFR", "FR", "CompoundFDP", "LinLog", "KK", "Circular", "Tree", "ISOM", "EfficientSugiyama" };
                graphLayout.LayoutAlgorithmType = this.ComboAlgorithm.SelectedItem.ToString();
                
                graphLayout.OverlapRemovalConstraint = AlgorithmConstraints.Must;   //.Automatic .Must .Skip
                graphLayout.OverlapRemovalAlgorithmType = "FSA";    /// FSA     OneWayFSA

                graphLayout.HighlightAlgorithmType = "Simple";      /// Simple  

                graphLayout.Relayout();
            }
        
            private void SetCustomLayout(int V, int E)
            {
                double density = (2.0000 * E) / (1.0000 * V * (V - 1.0000));              
                float mainFactor = (float) (Math.Log(V, 2) * density);

                /// Set distance between vertices
                int gapFactor = 50;
                float vertexGap = mainFactor * gapFactor;
                var overlapRemoval = new GraphSharp.Algorithms.OverlapRemoval.OverlapRemovalParameters();
                overlapRemoval.HorizontalGap = vertexGap;
                overlapRemoval.VerticalGap = vertexGap;
                this.graphLayout.OverlapRemovalParameters = overlapRemoval;
                
                /// Set vertices
                int vsizeFactor = 10;
                double esizeFactor = (V<10)?2:(V<20?1.5:(V<50?1:(V<100?0.6:0.3)));
                float vertexSize = mainFactor * vsizeFactor;
                double fsizeFactor = vsizeFactor * 0.8;
                foreach (Vertex  v in GraphToVisualize.Vertices)
                {
                    VertexControl VC = graphLayout.GetVertexControl(v);
                    VC.FontSize = vertexSize;
                    VC.FontWeight = FontWeights.UltraBold;
                    VC.MinWidth = vertexSize;
                    VC.MinHeight = vertexSize;
                    //VC.Padding = new Thickness(5, 10, 5, 10);
                }

                /// Set edges
                foreach (WeightedEdge<object> e in GraphToVisualize.Edges)
                {
                    EdgeControl EC = graphLayout.GetEdgeControl(e);
                    EC.FontSize = fsizeFactor;
                    EC.FontWeight = FontWeights.SemiBold;
                    EC.StrokeThickness = esizeFactor;
                }
            }

            private void ShowCyclesList()
            {
                /// Create List of cycles
                GraphCycles = new List<GraphCycle>();
                foreach (var item in GraphFromLibrary.minimumCycleBasis.Select((value, idx) => new { value, idx }))
                {
                    var cycle = item.value;
                    var index = item.idx;

                    double _weight = Math.Round(GraphFromLibrary.minimumCycleBasis_w[index], 4);
                    int _edges = cycle.Where(c => c == 1).ToArray().Length;
                    StringBuilder _cycle = new System.Text.StringBuilder();

                    foreach (var edg in cycle.Select((value, idx) => new { value, idx }))
                    {
                        if(edg.value == 1)
                        {
                            string v1 = GraphFromLibrary.EdgesReference[edg.idx].Item1.ToString();
                            string v2 = GraphFromLibrary.EdgesReference[edg.idx].Item2.ToString();
                            _cycle.Append("(" + v1 + "," + v2 + ") ");
                        }
                    }

                    GraphCycle gc = new GraphCycle { Weight = _weight, Edges = _edges, Cycle = _cycle.ToString()};
                    GraphCycles.Add(gc);
                }

                /// Assign List to DataGrid
                this.DataGrid1.ItemsSource = GraphCycles;
                this.DataGrid1.Columns.Clear();
                this.DataGrid1.AutoGenerateColumns = false;  
                this.DataGrid1.AutoGenerateColumns = true;  
                this.DataGrid1.Columns[0].Header = "W";
                this.DataGrid1.Columns[1].Header = "E";
                this.DataGrid1.Columns[2].Header = "Cycle";
                this.DataGrid1.Items.Refresh();

                /// Autosize columns
                foreach (DataGridColumn column in DataGrid1.Columns)
                {
                    column.Width = new DataGridLength(1.0, DataGridLengthUnitType.SizeToCells);
                }
            }

            private void ShowCliquesList()
            {
                /// Create List of cliques
                GraphCliques = new List<GraphClique>();
                ///foreach (var item in GraphFromLibrary.maximalCliqueSet.Select((value, idx) => new { value, idx }))
                for (int idx = GraphFromLibrary.maximalCliqueSet.Length; idx --> 0;)
                {
                    var clique = GraphFromLibrary.maximalCliqueSet[idx];;///item.value;
                    var index = idx;///item.idx;

                    int[] sortedClique = clique.OrderBy(i => i).ToArray();
                    int _number = sortedClique.Length;
                    string _vertices = string.Join(", ", sortedClique);
                
                    GraphClique gc = new GraphClique { Number = _number, Vertices = _vertices};
                    GraphCliques.Add(gc);
                }

                /// Assign List to DataGrid

                this.DataGrid1.ItemsSource = GraphCliques;
                this.DataGrid1.Columns.Clear();
                this.DataGrid1.AutoGenerateColumns = false;  
                this.DataGrid1.AutoGenerateColumns = true;                
                this.DataGrid1.Columns[0].Header = "K";
                this.DataGrid1.Columns[1].Header = "Vertices";
                this.DataGrid1.Items.Refresh();

                /// Autosize columns
                foreach (DataGridColumn column in DataGrid1.Columns)
                {
                    column.Width = new DataGridLength(1.0, DataGridLengthUnitType.SizeToCells);
                }
            }
            
            /*public void ChangeEdgeColor(object src, object dst, Color newColor)
            {
                WeightedEdge<object> edgeA = GraphToVisualize.Edges.ToList().Find(
                                    x => x.Source.ToString() == src.ToString() && x.Target.ToString() == dst.ToString());
                
                WeightedEdge<object> edgeB = GraphToVisualize.Edges.ToList().Find(
                                    x => x.Source.ToString() == dst.ToString() && x.Target.ToString() == src.ToString());

                if (edgeA is WeightedEdge<object>)
                {
                    edgeA.Color = newColor;
                    Application.Current.Dispatcher.Invoke(new Action(() => { GraphToVisualize.RemoveEdge(edgeA); }));
                    Application.Current.Dispatcher.Invoke(new Action(() => { GraphToVisualize.AddEdge(edgeA); }));
                }

                if (edgeB is WeightedEdge<object>)
                {
                    edgeB.Color = newColor;
                    Application.Current.Dispatcher.Invoke(new Action(() => { GraphToVisualize.RemoveEdge(edgeB); }));
                    Application.Current.Dispatcher.Invoke(new Action(() => { GraphToVisualize.AddEdge(edgeB); }));
                }
            }*/
            
            /*public void ChangeVertexColor(object v, Color newColor)
            {
                Vertex vertex = GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == v.ToString()) as Vertex;

                if (vertex is Vertex)
                {
                    vertex.Color = newColor;

                    var edges = new List<WeightedEdge<object>>(GraphToVisualize.Edges.Where( e =>
                        ((Vertex)e.Source == vertex) || ((Vertex) e.Target == vertex)));

                    Application.Current.Dispatcher.Invoke(new Action(() => { GraphToVisualize.RemoveVertex(vertex); }));
                    Application.Current.Dispatcher.Invoke(new Action(() => { GraphToVisualize.AddVertex(vertex); }));

                    foreach (var edge in edges)
                    {
                        if (edge is WeightedEdge<object>)
                        {
                            Application.Current.Dispatcher.Invoke(new Action(() => { GraphToVisualize.AddEdge(edge); }));
                        }
                    }
                }
            }*/

            #endregion

            #region EVENT HANLDERS

            /// On windows resize
            private void root_SizeChanged(object sender, SizeChangedEventArgs e)
            {
                var prevSize = e.PreviousSize;
                var newSize = e.NewSize;
                var resizeFactor = e.NewSize.Width / prevSize.Width;
                if(e.NewSize.Width > 0 && prevSize.Width > 0)
                {
                    DataGrid1.Width *= resizeFactor;
                    
                    Thickness m = BtnShowAll.Margin;
                    m.Left = DataGrid1.Width + 10;
                    BtnShowAll.Margin = m;

                    BorderFilters.Width *= resizeFactor;
                    SliderFilters.Width = (SliderFilters.Width * resizeFactor) - 30;
                }
            }

            /// On windows exit
            private void root_Closing(object sender, CancelEventArgs e)
            {
                //if(graphLayout.Dispatcher.Thread.ThreadState == ThreadState.Running)
                graphLayout.Dispatcher.Thread.Interrupt();
            }

            /// DataGrid selection changed
            private void DataGrid1_selectedCellsChanged(object sender, SelectedCellsChangedEventArgs  e)
            {
                /// Using selected cell
                //var item = e.AddedCells[0];
                //var col = item.Column as DataGridColumn;
                //var fc = col.GetCellContent(item.Item);
                
                var selIdx = DataGrid1.SelectedIndex;
                
                if(selIdx >= 0)
                {
                    int[] data;
                    int type, size, idx;
                    if (ShowingNow == 1)
                    {
                        type = 0;
                        idx = selIdx;                                                   /// because shown in ASCcending order
                        data = GraphFromLibrary.minimumCycleBasis[idx];
                        size = data.Length;
                    }
                    else
                    {
                        type = 2;
                        idx = (GraphFromLibrary.maximalCliqueSet.Length-1) - selIdx;    /// because shown in DESCending order
                        data = GraphFromLibrary.maximalCliqueSet[idx];
                        size = data.Length;
                        }

                    /// Run normally
                    MainWindow.HighlightPath(ref data, type, size, 1, false, true, Colors.Transparent);

                    /// Run in a new thread, with max priority (slows down the timed-out visualization)
                    /*Task.Run(()=> { Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Send,
                        new Action(() =>
                        {
                            HighlightPath(ref data, type, size, 1);

                        }));
                    });*/
                }
            }

            /// DataGridRow selected
            private void DataGrid1_RowClicked(object sender, RoutedEventArgs e)
            {
                //var index = (sender as DataGridRow).GetIndex();
                ///Start highlight
                Console.WriteLine("Mouse Down on Grid");
                MouseDownOnGrid = true;
            }

            /// DataGrid selected
            private void DataGrid1_Selected(object sender, RoutedEventArgs e)
            {
                /// Restore scroll position to left
                Border border = VisualTreeHelper.GetChild(DataGrid1, 0) as Border;
                ScrollViewer scrollViewer = border.Child as ScrollViewer;
                scrollViewer.ScrollToLeftEnd();

                /// Restore scroll position to left
                scrollViewer.ScrollToLeftEnd();
            }
            
            /// Mouse up on grid
            private void DataGrid1_MouseUp(object sender, MouseEventArgs e)
            {
                ///End highlight
                Console.WriteLine("Mouse Up on Grid");
                MouseDownOnGrid = false;
            }

            /// Slide across filtrations
            private void SliderFilters_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
            {
                if (SliderFilters != null)
                {
                    int filterN = Convert.ToInt32(SliderFilters.Value);
                    if(filterN >= 1)
                        FiltrationLabel.Content = MainWindow.AnimatedFiltration(0, filterN);
                }
            }

            #endregion
    }

    #region OTHER-CLASSES

    #endregion
}