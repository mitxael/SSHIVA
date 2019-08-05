/// ***********************************************************
///	File	: MainWindow.xaml.css
///	About	: Main GUI
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using System.Timers;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;
using CHoleR.Helper;
using GRAPH;
using QuickGraph;
using GraphSharp;
using GraphSharp.Controls;

namespace CHoleR
{
    public partial class MainWindow : Window
    {
        #region CLASS VARIABLES

        public static GraphManager GraphFromLibrary;
        public static BidirectionalWeightedGraph<object, WeightedEdge<object>> GraphToVisualize;
        public static GraphLayout graphLayout;
        public static List<TaggedVertex> HighlightedVertices = new List<TaggedVertex>();
        public static List<TaggedEdge> HighlightedEdges = new List<TaggedEdge>();
        private readonly BackgroundWorker backgroundWorker1 = new BackgroundWorker();
        public Boolean pendingChanges;
        
        #endregion

        #region WINDOW STARTUP

        public MainWindow()
        {
            InitializeComponent();

            InitializeBackgroundWorker();

            InitializeGUI();

            Resources["AppTitle"] = "Graph Algorithms";
            Resources["AppDirectory"] = System.IO.Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location);
            Resources["InputDirectory"] = Application.Current.MainWindow.Resources["AppDirectory"].ToString() + "\\input\\";
            Resources["OutputDirectory"] = Application.Current.MainWindow.Resources["AppDirectory"].ToString() + "\\output\\";
        }

        private void InitializeGUI()
        {
            this.RadioBSource1.IsChecked = true;

            this.RadioBAlgo3.IsChecked = true;

            this.ComboVerbosity.Items.Add("Minimum (faster)");
            this.ComboVerbosity.Items.Add("Normal");
            this.ComboVerbosity.Items.Add("Maximum (slower)");
            this.ComboVerbosity.SelectedIndex = 0;
        }

        #endregion

        #region BUTTON ACTIONS

        /// BTN EXECUTE ALGORITHM
        private void BtnAlgorithm_Click(object sender, RoutedEventArgs e)
        {
            /// Disable button to prevent multi-execution
            Btn_Start.IsEnabled = false;

            /// Reinitialize
            this.TextBox1.Clear();
            ProgressBar1.Value = -1;
            TextBlock1.Text = 0 + "%";

            if(pendingChanges == true)
                buildGraph(false);

            /// Console.SetOut(TextWriter.Synchronized(new TextBoxWriter(TextBox1)));
            try
            {
                ExecuteAlgorithm_Tasks();
            }
            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
            }
        }

        /// BTN SHOW GRAPH
        private void BtnViewGraph_Click(object sender, RoutedEventArgs e)
        {
            if(pendingChanges == true)
                buildGraph(false);

            /// Create the graph visualization
            CreateGraphToVisualize(0);

            var graphView = new GraphView(ref GraphFromLibrary, ref GraphToVisualize);

            try
            {
                /// SHOW() WITH "SEND" (...<Background<Input<Loaded<Render<DataBind<Normal<Send) PRIORITY              
                Task.Run(()=> { Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Send,
                        new Action(() =>
                        {
                            try
                            {
                                ///graphView.Show();
                                graphView.ShowDialog();
                            }
                            catch (Exception exception)
                            {
                                Console.WriteLine(exception);
                                ///throw;
                            } 
                        }));
                });
            }
            catch (ThreadInterruptedException ex)
            {
                Console.WriteLine("Thread interrupted! (Main) - " + ex);
            }
        }

        /// BTN SAVE GRAPH
        private void BtnSaveGraph_Click(object sender, RoutedEventArgs e)
        {
            MessageBoxResult messageBoxResult = System.Windows.MessageBox.Show("\nDo you confirm saving the graph?\n",
                                                                                "Save", System.Windows.MessageBoxButton.YesNo, MessageBoxImage.Question);
            if (messageBoxResult == MessageBoxResult.Yes)
            {
                try
                {
                    /// Save graph
                    if(pendingChanges == false)
                        GraphFromLibrary.saveGraph(GraphFromLibrary.source, Convert.ToSByte('A'));
                    else
                        buildGraph(true);
                }
                catch (Exception ex)
                {
                    Console.WriteLine(ex.Message);
                }
                finally
                {
                    pendingChanges = false;
                    MessageBoxResult msg = MessageBox.Show("Graph saved into \\input\\ folder as: " + GraphFromLibrary.source, "Save completed", MessageBoxButton.OK, MessageBoxImage.Warning);
                }
            }
        }

        /// BTN BENCHMARK ALGORITHMS
        private async void BtnBenchmark_Click(object sender, RoutedEventArgs e)
        {
            ///TODO: Select graphs to benchmark (new window w/listOfCheckboxes) (also for algorithm)
            
            String GraphList = "";
            foreach (var item in ComboGraph.Items)
            {
                GraphList = String.Join("\n", GraphList, item.ToString());
            }

            string algoType = (RadioBAlgo1.IsChecked == true) ? "Minimal Cycle Basis" : "Maximal Cliques";

            MessageBoxResult messageBoxResult = System.Windows.MessageBox.Show("\nThis operation may take several minutes." +
                                                                               " Do you confirm benchmarking all files in output folder?\n" +
                                                                               GraphList + "\n",
                                                                               "Benchmark for "+algoType, System.Windows.MessageBoxButton.YesNo, MessageBoxImage.Question);
            if (messageBoxResult == MessageBoxResult.Yes)
            {
                /// Disable button to prevent multi-execution
                Btn_Benchmark.IsEnabled = false;

                /// Show waiting bar
                TextBlock_Wait.Visibility = Visibility.Visible;
                ForceRefresh(TextBlock_Wait);
                

                /// Reinitialize
                this.TextBox1.Clear();
                ProgressBar1.Value = -1;
                TextBlock1.Text = 0 + "%";

                /// Activate the list of files
                RadioBSource1.IsChecked = true;

                /// Set verbosity to minimum
                ComboVerbosity.SelectedIndex = 0;

                ///Create the table of results
                var results = new List<Tuple<string, string, double, double>>();

                /// For each graph file
                int i_max = ComboGraph.Items.Count;
                //int i_max = 1;
                for (int i = 0; i < i_max; ++i)
                {
                    ComboGraph.SelectedIndex = i;

                    buildGraph(false);

                    /// Perform each algorithm
                    int j_max = ComboAlgorithm.Items.Count;
                    //int j_max = 4; /// {Horton, DePina, Hybrid, Amaldi}
                    double[] benchmark = new double[j_max];
                    var tasks = new List<Task>();             

                    for (int j = 0; j < j_max; ++j)
                    {
                        ComboAlgorithm.SelectedIndex = j;

                        try
                        {                          
                            /// Serial execution of algorithms for any graph
                            //benchmark[j] = await ExecuteAlgorithm_Tasks();

                            /// Parallel execution of "n" algorithms for the same graph
                            var task = ExecuteAlgorithm_Tasks();
                            benchmark[j] = await task;
                            tasks.Add(task);
                        }
                        catch (Exception ex)
                        {
                            MessageBox.Show(ex.Message);
                        }
                        finally
                        ///Create a table and insert: G1.source, algorithmName, algorithmTime
                        {
                            /// Display graph's name only for the first line of the block
                            string label = (j==0)? GraphFromLibrary.source : "";  
  
                            /// Set RESULTS, including either MCB weight or CLIQUES number
                            if(RadioBAlgo1.IsChecked == true)
                                results.Add(new Tuple<string, string, double, double>(label, ComboAlgorithm.SelectedItem.ToString(), GraphFromLibrary.minimumCycleBasis_w.Sum(), benchmark[j]));
                            else
                                results.Add(new Tuple<string, string, double, double>(label, ComboAlgorithm.SelectedItem.ToString(), GraphFromLibrary.maximalCliqueSet.Length, benchmark[j]));

                            /// Interval to prevent contemporary access to the output files
                            Thread.Sleep(100);
                        }
                    }

                    /// Wait for all algorithms to complete before processing next graph
                    await Task.WhenAll(tasks);
                }

                /// Remove waiting bar
                TextBlock_Wait.Visibility = Visibility.Hidden;

                /// Re-enable button
                Btn_Benchmark.IsEnabled = true;

                /// At the end, save results to a .csv file
                MessageBoxResult msg = MessageBox.Show("Results saved into \\output\\ folder as benchmark.csv", "Benchmark completed", MessageBoxButton.OK, MessageBoxImage.Warning);
                SaveToExcel(results, true);
            }
        }

        #endregion

        #region ALGORITHM EXECUTION

        private async Task<double> ExecuteAlgorithm_Tasks()
        {
            /// Result
            double algorithmTime = -1;

            /// Build Graph
            //buildGraph();
            
            /// Execute algorithm
            try
            {
                /// Execution parameters
                int algoType = (RadioBAlgo1.IsChecked == true) ? 0 : ((RadioBAlgo2.IsChecked == true)? 4 : 8);
                int algorithm = int.Parse(this.ComboAlgorithm.SelectedIndex.ToString()) + 1 + algoType;
                int verbosity = int.Parse(this.ComboVerbosity.SelectedIndex.ToString());
                int foretime = (GraphFromLibrary.E * GraphFromLibrary.V / 10000) * 2;
                bool hiresClock = true;
                bool saveToFile = true;

                /// Algorithm parameters
                List<Tuple<string, string>> algoParams = new List<Tuple<string, string>>();
                algoParams.Add(new Tuple<string, string>(LabelAlgorithm.Content.ToString(), ComboAlgorithm.SelectedIndex.ToString()));
                algoParams.Add(new Tuple<string, string>(LabelParam1.Content.ToString(), ComboParam1.SelectedIndex.ToString()));
                algoParams.Add(new Tuple<string, string>(LabelParam2.Content.ToString(), ComboParam2.SelectedIndex.ToString()));
                algoParams.Add(new Tuple<string, string>(LabelParam3.Content.ToString(), ComboParam3.SelectedIndex.ToString()));
                algoParams.Add(new Tuple<string, string>(LabelParam4.Content.ToString(), ComboParam4.SelectedIndex.ToString()));

                switch (algorithm)
                {
                    /// Horton
                    case 1:
                        foretime /= 1;
                        break;
                    /// Depina
                    case 2:
                        foretime /= 2;
                        break;
                    /// Hybrid
                    case 3:
                        foretime /= 3;
                        break;
                    /// Amaldi
                    case 4:
                        foretime /= 4;
                        break;
                    /// Criccaldi
                    case 5:
                        foretime /= 4;
                        GraphFromLibrary.K = 3;
                        break;
                    /// BronKerbosch v1
                    case 6:
                        foretime /= 6;
                        GraphFromLibrary.K = 3;
                        break;
                    /// BronKerbosch v2
                    case 7:
                        foretime /= 6;
                        GraphFromLibrary.K = 3;
                        break;
                    /// BronKerbosch v3
                    case 8:
                        foretime /= 6;
                        GraphFromLibrary.K = 3;
                        break;
                    /// PH Javaplex, Gudhi, Phat, Dipha, Perseus
                    case 9: case 10: case 11: case 12: case 13:
                        ///TODO: According to density
                        foretime /= 1;
                        break;
                    default:
                        foretime /= 1;
                        break;
                }
                UtilsManager start = new UtilsManager();
                
                try
                {
                    ///Show debug-info
                    /*MessageBoxResult messageBoxResult = System.Windows.MessageBox.Show("algorithm: " + algorithm,
                                    "Debug info", System.Windows.MessageBoxButton.OK, MessageBoxImage.Hand);*/
                    
                    ///Run tasks
                    Task t1 = Task.Run(() => ShowProgress(foretime, ref algorithmTime));
                    Task t2 = Task.Run(() => algorithmTime = start.measurePerformance(ref GraphFromLibrary, algorithm, verbosity, hiresClock, saveToFile, algoParams));

                    ///TODO: Add timeout for execution
                   /* int timeout = 1*1000; /// timeout in milliseconds
                    var linkedTokenSource = CancellationTokenSource.CreateLinkedTokenSource(new CancellationTokenSource(timeout).Token);

                    Task t2 = Task.Run(() =>
                    {
                        algorithmTime = start.measurePerformance(ref G1, algorithm, verbosity, 0, 1);
                        linkedTokenSource.Token.ThrowIfCancellationRequested();
                    }, linkedTokenSource.Token);
 
                    if ( await Task.WhenAny(t1, Task.Delay(timeout)) != t2)
                    {
                        linkedTokenSource.Cancel();
                        TextBox1.Text = "Timeout! The limit of " + timeout + "ms have been reached.";
                        Thread.CurrentThread.Interrupt();
                        t2.Dispose();
                        throw new TimeoutException(TextBox1.Text);
                    }*/

                    ///Execute task with highest priority
                    /*Task t2 = Task.Run(()=> { Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Send,
                        new Action(() =>
                        {
                            algorithmTime = start.measurePerformance(ref G1, algorithm, verbosity, 0, 1);

                        }));
                    });*/
                    
                    await Task.WhenAll(t1, t2);
                }
                catch (Exception ex)
                {
                    //MessageBox.Show(ex.Message);
                    Console.WriteLine(ex.Message);
                    throw;
                }

            }
            finally
            {
                try
                {
                    /// release all resources
                    //start.terminateExecution();

                    if(algorithmTime > 0)
                    {
                        /// Determine the last "result_"n".txt" file to read
                        string outputDirectory = Application.Current.MainWindow.Resources[$"OutputDirectory"].ToString();
                        var here = Directory.EnumerateFiles(outputDirectory, "result_*.txt");
                        int counter = 0;
                        foreach (string file in here)
                        {
                            string fileName = System.IO.Path.GetFileName(file); //"result_xyz.txt"
                            var filenameChunks = fileName.Split('_','.'); // one or more separator characters
                            var filenameNumChunks = filenameChunks.Where(chunk => chunk.All(char.IsDigit));
                            var lastNumber = Convert.ToInt32(String.Concat(filenameNumChunks));
                            if (lastNumber > counter) counter = lastNumber;
                        }
                        String result_filename = "result_"+ (counter).ToString() + ".txt";

                        /// Get result into TextBox
                        String result = outputDirectory + result_filename;
                        TextBox1.Text = File.ReadAllText(result);
                    }
                }
                catch (Exception ex)
                {
                    MessageBox.Show(ex.Message);
                }
                //TextBox1.AppendText($"\n\n **** CPU TIME: {algorithmTime:0.000000##} seconds ****");
                TextBox1.ScrollToEnd();

                /// Enable "Start" button to allow another measurement
                Btn_Start.IsEnabled = true;
            }

            return algorithmTime;
        }

        /*private void ExecuteAlgorithm_Threads()
        {

            //Prepare Graph
            int algorithm = int.Parse(this.ComboAlgorithm.SelectedIndex.ToString()) + 1;
            int verbosity = int.Parse(this.ComboVerbosity.SelectedIndex.ToString());
            string filename = this.ComboGraph.SelectedItem.ToString();
            GraphManager G1 = new GraphManager(0);
            G1.buildGraph('A', filename);

            //Run Tasks!
            ParallelOptions po = new ParallelOptions();
            po.MaxDegreeOfParallelism = 10;
            Parallel.Invoke(() =>
               {
                   //Show Progress
                   int V = Int32.Parse(File.ReadLines("G:\\DOCUMENTS\\DEVELOPMENT\\C++\\graph-algorithms\\resources\\" + filename).First());
                   int E = Int32.Parse(File.ReadLines("G:\\DOCUMENTS\\DEVELOPMENT\\C++\\graph-algorithms\\resources\\" + filename).Skip(1).FirstOrDefault());
                   int foretime = E*V/1000;
                   var t1 = new Thread(() => ShowProgress(foretime));    //lamba usage to pass parameters
                   t1.Start();
               },

               () =>
               {
                   //MeasureAlgorithm
                   try
                   {
                       Console.WriteLine("+++++ t2 started");
                       UtilsManager start = new UtilsManager();
                       AlgorithmTime = start.measurePerformance(G1, algorithm, verbosity, 0);
                       TextBox1.Dispatcher.BeginInvoke(System.Windows.Threading.DispatcherPriority.Normal,
                           new Action(delegate ()
                           {
                               TextBox1.AppendText($"{AlgorithmTime:0.000000##}");
                           })
                       );
                   }
                   catch (Exception Ex)
                   {
                       this.TextBox1.AppendText("Error in system: " + Ex.Message + " ---> " + Ex.StackTrace);
                   }
                   Console.WriteLine("+++++ t2 ended");
               }
           );
        }*/

        private void ShowProgress(int delay, ref double algorithmTime)
        {
            for (int n = 0; n < 100 && algorithmTime < 0; ++n)
            {
                Thread.Sleep(delay);    /// Blocks the thread (from doing anything) for "delay" milliseconds
                ProgressBar1.Dispatcher.Invoke(System.Windows.Threading.DispatcherPriority.Send,
                    new Action(delegate()
                        {
                            ProgressBar1.Value = n;
                        }
                    ));

                TextBlock1.Dispatcher.Invoke(System.Windows.Threading.DispatcherPriority.Send,
                    new Action(delegate()
                        {
                            TextBlock1.Text = n + "%";
                        }
                    ));
                if (n == 90 && algorithmTime <= 0)
                {
                    n -= 1;
                }
            }

            //Complete in case of early exit
            ProgressBar1.Dispatcher.Invoke(System.Windows.Threading.DispatcherPriority.Send,
                new Action(delegate ()
                    {
                        ProgressBar1.Value = 100;
                    }
                ));
            TextBlock1.Dispatcher.Invoke(System.Windows.Threading.DispatcherPriority.Send,
                new Action(delegate ()
                    {
                        TextBlock1.Text = "100%";
                    }
                ));
        }

        private void SetParamsOnGUI()
        {
            int algoType = (RadioBAlgo1.IsChecked == true) ? 0 : ((RadioBAlgo2.IsChecked == true)? 4 : 8);
            int algorithm = int.Parse(this.ComboAlgorithm.SelectedIndex.ToString()) + 1 + algoType;

            switch (algorithm)
            {
                /// Horton
                case 1:
                {
                    LabelAlgorithm.Content = "Algorithm";
                    
                    LabelParam1.Visibility = Visibility.Hidden;
                    ComboParam1.Visibility = Visibility.Hidden;
                    LabelParam2.Visibility = Visibility.Hidden;
                    ComboParam2.Visibility = Visibility.Hidden;
                    LabelParam3.Visibility = Visibility.Hidden;
                    ComboParam3.Visibility = Visibility.Hidden;
                    LabelParam4.Visibility = Visibility.Hidden;
                    ComboParam4.Visibility = Visibility.Hidden;

                    break;
                }

                /// Depina
                case 2:
                {
                    LabelAlgorithm.Content = "Algorithm";

                    LabelParam1.Content = "Cycles finder";
                    ComboParam1.Items.Clear();
                    ComboParam1.Items.Add("SignedGraph");
                    ComboParam1.Items.Add("HortonSpace");
                    ComboParam1.SelectedIndex = 0;

                    LabelParam2.Visibility = Visibility.Hidden;
                    ComboParam2.Visibility = Visibility.Hidden;
                    LabelParam3.Visibility = Visibility.Hidden;
                    ComboParam3.Visibility = Visibility.Hidden;
                    LabelParam4.Visibility = Visibility.Hidden;
                    ComboParam4.Visibility = Visibility.Hidden;

                    break;
                }

                /// Hybrid
                case 3:
                {
                    LabelAlgorithm.Content = "Algorithm";

                    LabelParam1.Visibility = Visibility.Hidden;
                    ComboParam1.Visibility = Visibility.Hidden;
                    LabelParam2.Visibility = Visibility.Hidden;
                    ComboParam2.Visibility = Visibility.Hidden;
                    LabelParam3.Visibility = Visibility.Hidden;
                    ComboParam3.Visibility = Visibility.Hidden;
                    LabelParam4.Visibility = Visibility.Hidden;
                    ComboParam4.Visibility = Visibility.Hidden;

                    break;
                }

                /// Amaldi
                case 4:
                {
                    LabelAlgorithm.Content = "Algorithm";

                    LabelParam1.Visibility = Visibility.Visible;
                    LabelParam1.Content = "Only isometric";
                    ComboParam1.Visibility = Visibility.Visible;
                    ComboParam1.Items.Clear();
                    ComboParam1.Items.Add("Yes");
                    ComboParam1.Items.Add("No");
                    ComboParam1.SelectedIndex = 0;

                    LabelParam2.Visibility = Visibility.Hidden;
                    ComboParam2.Visibility = Visibility.Hidden;
                    LabelParam3.Visibility = Visibility.Hidden;
                    ComboParam3.Visibility = Visibility.Hidden;
                    LabelParam4.Visibility = Visibility.Hidden;
                    ComboParam4.Visibility = Visibility.Hidden;

                    break;
                }

                /// Criccaldi
                case 5:
                {
                    LabelAlgorithm.Content = "Algorithm";

                    LabelParam1.Visibility = Visibility.Visible;
                    LabelParam1.Content = "Minimum size (k)";
                    ComboParam1.Visibility = Visibility.Visible;
                    ComboParam1.Items.Clear();
                    int graphOrder = 3;
                    for(int i = 0; i < graphOrder; ++i)
                        ComboParam1.Items.Add((i+1).ToString());
                    ComboParam1.SelectedIndex = 2;

                    LabelParam2.Visibility = Visibility.Hidden;
                    ComboParam2.Visibility = Visibility.Hidden;
                    LabelParam3.Visibility = Visibility.Hidden;
                    ComboParam3.Visibility = Visibility.Hidden;
                    LabelParam4.Visibility = Visibility.Hidden;
                    ComboParam4.Visibility = Visibility.Hidden;

                    break;
                }

                /// BronKerbosch v1
                case 6:
                {
                    LabelAlgorithm.Content = "Algorithm";

                    LabelParam1.Visibility = Visibility.Visible;
                    LabelParam1.Content = "Minimum size (k)";
                    ComboParam1.Visibility = Visibility.Visible;
                    ComboParam1.Items.Clear();
                    int graphOrder = 3;
                    for(int i = 0; i < graphOrder; ++i)
                        ComboParam1.Items.Add((i+1).ToString());
                    ComboParam1.SelectedIndex = 2;

                    LabelParam2.Visibility = Visibility.Hidden;
                    ComboParam2.Visibility = Visibility.Hidden;
                    LabelParam3.Visibility = Visibility.Hidden;
                    ComboParam3.Visibility = Visibility.Hidden;
                    LabelParam4.Visibility = Visibility.Hidden;
                    ComboParam4.Visibility = Visibility.Hidden;

                    break;
                }

                /// BronKerbosch v2
                case 7:
                {
                    LabelAlgorithm.Content = "Algorithm";

                    LabelParam1.Visibility = Visibility.Visible;
                    LabelParam1.Content = "Minimum size (k)";
                    ComboParam1.Visibility = Visibility.Visible;
                    ComboParam1.Items.Clear();
                    int graphOrder = 3;
                    for(int i = 0; i < graphOrder; ++i)
                        ComboParam1.Items.Add((i+1).ToString());
                    ComboParam1.SelectedIndex = 2;

                    LabelParam2.Visibility = Visibility.Hidden;
                    ComboParam2.Visibility = Visibility.Hidden;
                    LabelParam3.Visibility = Visibility.Hidden;
                    ComboParam3.Visibility = Visibility.Hidden;
                    LabelParam4.Visibility = Visibility.Hidden;
                    ComboParam4.Visibility = Visibility.Hidden;

                    break;
                }

                /// BronKerbosch v3
                case 8:
                {
                    LabelAlgorithm.Content = "Algorithm";

                    LabelParam1.Visibility = Visibility.Visible;
                    LabelParam1.Content = "Minimum size (k)";
                    ComboParam1.Visibility = Visibility.Visible;
                    ComboParam1.Items.Clear();
                    int graphOrder = 3;
                    for(int i = 0; i < graphOrder; ++i)
                        ComboParam1.Items.Add((i+1).ToString());
                    ComboParam1.SelectedIndex = 2;

                    LabelParam2.Visibility = Visibility.Hidden;
                    ComboParam2.Visibility = Visibility.Hidden;
                    LabelParam3.Visibility = Visibility.Hidden;
                    ComboParam3.Visibility = Visibility.Hidden;
                    LabelParam4.Visibility = Visibility.Hidden;
                    ComboParam4.Visibility = Visibility.Hidden;

                    break;
                }

                /// PH Javaplex, Gudhi, PHAT, DiPHA, Perseus
                case 9: case 10: case 11: case 12:
                {
                    LabelAlgorithm.Content = "PH engine";
                    ComboAlgorithm.Visibility = Visibility.Visible;
                    //ComboAlgorithm...

                    LabelParam1.Visibility = Visibility.Visible;
                    LabelParam1.Content = "Filtration";
                    ComboParam1.Visibility = Visibility.Visible;
                    ComboParam1.Items.Clear();
                    ComboParam1.Items.Add("AllCliques");
                    if(algorithm == 9)
                    {
                        ComboParam1.Items.Add("Graph");
                        ComboParam1.Items.Add("MaxCliques");
                        ComboParam1.Items.Add("MinCycles");
                    }
                    ComboParam1.SelectedIndex = 0;

                    LabelParam2.Visibility = Visibility.Visible;
                    LabelParam2.Content = "Intervals";
                    ComboParam2.Visibility = Visibility.Visible;
                    ComboParam2.Items.Clear();
                    ComboParam2.Items.Add("Only infinite");
                    ComboParam2.Items.Add("All intervals");
                    ComboParam2.SelectedIndex = 0;

                    LabelParam3.Visibility = Visibility.Visible;
                    LabelParam3.Content = "Dimensions";
                    ComboParam3.Visibility = Visibility.Visible;
                    ComboParam3.Items.Clear();
                    ComboParam3.Items.Add("Show all");
                    
                    ComboParam3.Items.Add("Show dim-0 only");
                    ComboParam3.Items.Add("Show dim-1 only");
                    ComboParam3.Items.Add("Show dim-2 only");
                    ComboParam3.SelectedIndex = 0;

                    if(algorithm == 9)
                    {
                        LabelParam4.Visibility = Visibility.Visible;
                        LabelParam4.Content = "Holes";
                        ComboParam4.Visibility = Visibility.Visible;
                        ComboParam4.Items.Clear();
                        ComboParam4.Items.Add("Save as graph cycles");
                        ComboParam4.Items.Add("Do not save");
                        ComboParam4.SelectedIndex = 0;
                    }
                    else if (algorithm == 10)
                    {                       
                        LabelParam4.Visibility = Visibility.Hidden;
                        ComboParam4.Visibility = Visibility.Hidden;
                    }

                    break;
                }

                /// Default
                default:
                {
                    break;
                }
            }
        }

        #endregion

        #region TOOL FUNCTIONS

        public class TextBoxWriter : TextWriter
        {
            TextBox _output = null;

            public TextBoxWriter(TextBox output)
            {
                _output = output;
            }

            public override void Write(char value)
            {
                base.Write(value);
                _output.AppendText(value.ToString());
            }

            public override Encoding Encoding
            {
                get { return System.Text.Encoding.UTF8; }
            }
        }

        private void NumberValidationTextBox(object sender, TextCompositionEventArgs e)
        {
            /// Integer (default)
            Regex regex = new Regex("[^0-9]{1,4}$");

            /// Float (density graph)
            if( RadioBSource4.IsChecked == true && ((TextBox)sender).Name == "Text_b")
                regex = new Regex("^[a-zA-Z]{1,3}$");
            //regex = new Regex("^[0-9]([.,][0-9]{1,3})?$");
            //regex = new Regex(@"^[0-9]*(?:\.[0-9]+)?$");

            e.Handled = regex.IsMatch(e.Text);
        }  

        public void SaveToExcel(List<Tuple<string, string, double, double>> results, Boolean openAfter)
        {
            string outputDirectory = Application.Current.MainWindow.Resources[$"OutputDirectory"].ToString();
            string path = outputDirectory + "\\benchmark.csv";
            using (var file = File.CreateText(path))
            {
                /// Wirte HEADER
                if(RadioBAlgo1.IsChecked == true)
                    file.WriteLine(String.Join("FILE", "ALGORITHM", "MCB_WEIGHT", "TIME (s)"));
                else
                    file.WriteLine(String.Join("FILE", "ALGORITHM", "CLIQUES_NUMBER", "TIME (s)"));

                foreach(var tuple in results)
                {
                    /// GraphSource, Algorithm, Time
                    //file.WriteLine(String.Join(", ", tuple.Item1, tuple.Item2, tuple.Item4));
                    /// GraphSource, Algorithm, MCB, Time
                    file.WriteLine(String.Join(", ", tuple.Item1, tuple.Item2, tuple.Item3, tuple.Item4));
                }
            }
            if(openAfter==true)
                Process.Start(path);
        }

        private delegate void NoArgDelegate();

        public static void ForceRefresh(DependencyObject obj)
        {
            try
            {
                obj.Dispatcher.Invoke(System.Windows.Threading.DispatcherPriority.ApplicationIdle,
                    (NoArgDelegate)delegate { });
            }
            catch (ThreadInterruptedException ex)
            {
                Console.WriteLine("Thread interrupted! (GraphLayout) - " + ex);
            }
            finally{}
        }
        
        /// Delay "x" seconds
        public static void Delay(int timeSec)
        {
            try
            {
                Thread.Sleep(timeSec);
            }
            catch (ThreadInterruptedException ex)
            {
                Console.WriteLine("Thread interrupted! (GraphLayout) - " + ex);
            }
        }

        ///  Get sequential colours
        public static Color GetSequentialColor(int index, bool random, int colorSelection)
        {
            Color seqColor;

            if(random)
            {
                Random rnd = new Random(); Byte[] b = new Byte[3]; rnd.NextBytes(b); 
                seqColor = Color.FromRgb(b[0],b[1],b[2]);
            }
            else
            {
                /// Prepare palette of colours
                Color[] colorsHard = { Colors.Red, Colors.Blue, Colors.Green, Colors.OrangeRed, Colors.Yellow, Colors.Purple };
                Color[] colorsMedium = { Colors.IndianRed, Colors.RoyalBlue, Colors.Olive, Colors.SaddleBrown, Colors.Goldenrod, Colors.HotPink };
                Color[] colorsSoft = { Colors.Tomato, Colors.Teal, Colors.LawnGreen, Colors.Salmon, Colors.Gold, Colors.MediumSlateBlue };
                Color[] colorsLight = { Colors.LightCoral, Colors.LightBlue, Colors.LightGreen, Colors.LightSalmon, Colors.LightGoldenrodYellow, Colors.LightPink };
                Color[] palette;
                
                /// Standard
                if(colorSelection == 0)
                    palette = colorsMedium.Union(colorsSoft).ToArray();
                /// Dark
                else if(colorSelection == 1)
                    palette = colorsHard.Union(colorsMedium).ToArray();
                /// Light
                else ///if(colorSelection == 2)
                    palette = colorsSoft.Union(colorsLight).ToArray();

                    
                seqColor = palette[(index + palette.Length) % palette.Length];
            }

            return seqColor;
        }

        public static Color HexToColor(string hex)
        {
            /*int r = Convert.ToInt32(hex.Substring(1,2),16);
            int g = Convert.ToInt32(hex.Substring(3,2),16);
            int b = Convert.ToInt32(hex.Substring(5,2),16);
            return Color.FromRgb((byte)(r & 0x000000FF), (byte)(g & 0x000000FF), (byte)(b & 0x000000FF));*/
            return (Color)ColorConverter.ConvertFromString(hex);
        }

        #endregion

        #region GRAPH METHODS
        
        private void buildGraph(Boolean save2file)
        {
            /// Initialize
            GraphFromLibrary = new GraphManager(0);

            /// File graph
            if (RadioBSource1.IsChecked == true)
            {
                string filename = this.ComboGraph.SelectedItem.ToString();
                GraphFromLibrary.buildGraph(filename);         
            }

            /// Hypercube graph
            if (RadioBSource2.IsChecked == true)
            {
                int v = this.ComboGraph.SelectedIndex + 2;
                GraphFromLibrary.buildGraph(save2file, v);
            }  
            
            /// Euclidean graph
            if (RadioBSource3.IsChecked == true)
            {
                int v = (int)SliderA.Value;
                int e = (int)SliderB.Value;
                GraphFromLibrary.buildGraph(save2file, v, e);
            }

            /// Density graph
            if (RadioBSource4.IsChecked == true)
            {
                int v = (int)SliderA.Value;
                double d = SliderB.Value;
                GraphFromLibrary.buildGraph(save2file, v, d);
            }

            /// After creation, no changes are pending
            pendingChanges = false;
        }

        public void CreateGraphToVisualize(int source)
        {
            /// INITIALIZE
            GraphToVisualize = new BidirectionalWeightedGraph<object, WeightedEdge<object>>();

            if (source == 0)
            {
                /// 0: READ FROM GRAPH
                if(GraphFromLibrary != null)
                {
                    for(int v = 0; v < GraphFromLibrary.V; ++v)
                    {
                        Color color = (Color) Application.Current.Resources["VertexNormal"];;
                        TaggedVertex vertex = new TaggedVertex(v, color);
                        GraphToVisualize.AddVertex(vertex);
                    }
                    foreach (var vertexList in GraphFromLibrary.adjList.Select((list, idx) => new { list, idx } ))
                    {
                        int src = vertexList.idx;
                        foreach (var edge in vertexList.list)
                        {
                            int dst = edge.Item1;
                            double weight = edge.Item2;
                            Color color = (Color) Application.Current.Resources["EdgeNormal"];;
                            
                            /// To improve visualization, add ONLY ONE bidirectional edge
                            if(src < dst)
                            {
                                ///Add colored edge (#7EADDB)
                                //WeightedEdge<object> wedge = new WeightedEdge<object>(src, dst , weight, color);
                                TaggedVertex v1 = (TaggedVertex) GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == src.ToString());
                                TaggedVertex v2 = (TaggedVertex) GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == dst.ToString());
                                TaggedEdge wedge = new TaggedEdge(v1, v2 , weight, color);
                                GraphToVisualize.AddEdge(wedge);
                            }
                        }
                    }
                }
                else
                {
                    MessageBoxResult msg = MessageBox.Show("There's no graph to plot.", "No graph was found", MessageBoxButton.OK, MessageBoxImage.Warning);
                }
            }
            else
            {
            /// 1: READ FROM FILE
                if (!String.IsNullOrEmpty(GraphFromLibrary.source))
                {
                    int V = 0;
                    int E = 0;

                    /// READ #VERTICES and #EDGES
                    int idx = 0;
                    string inputDirectory = Application.Current.MainWindow.Resources[$"InputDirectory"].ToString();
                    foreach (string line in File.ReadLines(inputDirectory + GraphFromLibrary.source).Take(2))
                    {
                        if (idx == 0)
                            V = Int32.Parse(line);
                        if (idx == 1)
                            E = Int32.Parse(line);
                        idx++;
                    }

                    /// ADD VERTICES
                    string[] vertices = new string[V];
                    for (int v = 0; v < V; ++v)
                    {
                        vertices[v] = v.ToString();
                        Color color = (Color) Application.Current.Resources["VertexNormal"];
                        TaggedVertex vertex = new TaggedVertex(v, color);
                        GraphToVisualize.AddVertex(vertex);
                    }

                    /// ADD EDGES AND WEIGHTS
                    var fileArray = File.ReadLines(inputDirectory + GraphFromLibrary.source).Skip(2).ToArray();
                    for (int e = 0; e < fileArray.Length; e++)
                    {
                        var columns = fileArray[e].Split(null); //default: space
                        int src = int.Parse(columns[0]);
                        int dst = int.Parse(columns[1]);
                        double weight = double.Parse(columns[2], System.Globalization.CultureInfo.InvariantCulture);
                        Color color = (Color) Application.Current.Resources["EdgeNormal"];

                        ///Add colored edge
                        TaggedVertex v1 = (TaggedVertex) GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == src.ToString());
                        TaggedVertex v2 = (TaggedVertex) GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == dst.ToString());
                        TaggedEdge wedge = new TaggedEdge(v1, v2 , weight, color);
                        GraphToVisualize.AddEdge(wedge);
                    }
                }
                else
                {
                    MessageBoxResult msg = MessageBox.Show("There's no graph to plot.", "No graph was found", MessageBoxButton.OK, MessageBoxImage.Warning);
                }
            }
        }

        private void ShowHideExtraFields(Boolean show)
        {
            if (show == true)
            {
                ComboGraph.Visibility = Visibility.Hidden;
                SliderA.Visibility = Visibility.Visible;
                Label_A.Visibility = Visibility.Visible;
                Label_Amin.Visibility = Visibility.Visible;
                Label_Amax.Visibility = Visibility.Visible;
                Label_Aval.Visibility = Visibility.Visible;
                SliderB.Visibility = Visibility.Visible;
                Label_B.Visibility = Visibility.Visible;
                Label_Bmin.Visibility = Visibility.Visible;
                Label_Bmax.Visibility = Visibility.Visible;
                Label_Bval.Visibility = Visibility.Visible;
                /// Euclidean graph
                if (RadioBSource3.IsChecked == true)
                {
                    //Label_Bval.Visibility = Visibility.Hidden;
                    Label_Bval.Content = "(edges)";

                    int V = Convert.ToInt32(SliderA.Value);
                    int maxE = V * (V - 1) / 2;
                    SliderB.Minimum = 0;
                    SliderB.Maximum = maxE;
                    SliderB.TickFrequency = 1.00;
                    SliderB.SmallChange = 1.00;
                    SliderB.LargeChange = SliderB.SmallChange;
                }
                /// Density graph
                if(RadioBSource4.IsChecked == true)
                {
                    //Label_Bval.Visibility = Visibility.Visible;
                    Label_Bval.Content = "(0 edges)";

                    SliderB.Minimum = 0.000;
                    SliderB.Maximum = 1.000;
                    SliderB.Value = SliderB.Minimum;
                    SliderB.TickFrequency = 0.001;
                    SliderB.SmallChange = 0.001;
                    SliderB.LargeChange = SliderB.SmallChange;
                }
            }
            else
            {
                ComboGraph.Visibility = Visibility.Visible;
                SliderA.Visibility = Visibility.Hidden;
                Label_A.Visibility = Visibility.Hidden;
                Label_Amin.Visibility = Visibility.Hidden;
                Label_Amax.Visibility = Visibility.Hidden;
                Label_Aval.Visibility = Visibility.Hidden;
                SliderB.Visibility = Visibility.Hidden;
                Label_B.Visibility = Visibility.Hidden;
                Label_Bmin.Visibility = Visibility.Hidden;
                Label_Bmax.Visibility = Visibility.Hidden;
                Label_Bval.Visibility = Visibility.Hidden; 
            }
        }

        public static void RemoveCustomHighlights()
        {
            if (MainWindow.HighlightedVertices != null && MainWindow.HighlightedVertices.Count > 0)
            {
                foreach (var v in HighlightedVertices)
                {
                    var vertex = GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == v.ToString()) as TaggedVertex;
                    vertex.Color = (Color)Application.Current.Resources["VertexNormal"];
                }
            }

            if (MainWindow.HighlightedEdges != null && MainWindow.HighlightedEdges.Count > 0)
            {
                foreach (var he in HighlightedEdges)
                {
                    he.Color = (Color)Application.Current.Resources["EdgeNormal"];

                    var v1 = GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == he.Source.ToString()) as TaggedVertex;
                    v1.Color = (Color)Application.Current.Resources["VertexNormal"];

                    var v2 = GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == he.Target.ToString()) as TaggedVertex;
                    v2.Color = (Color)Application.Current.Resources["VertexNormal"];
                }
            }
        }

        public static void HighlightVertex(string source, bool defaultColors, Color customColor)
        {
            ///My.Highlight
            Color VertexHighlightedColor;
            Color EdgeHighlightedColor;
            if (defaultColors)
            {
                VertexHighlightedColor = (Color) Application.Current.Resources["VertexHighlighted"];
                EdgeHighlightedColor = (Color) Application.Current.Resources["EdgeHighlighted"];
            }
            else
            {
                ///VertexHighlightedColor = EdgeHighlightedColor = GetSequentialColor(1, false);
                VertexHighlightedColor = EdgeHighlightedColor = customColor;
            }

            var vertex = GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == source) as TaggedVertex;
            if (vertex.Color != VertexHighlightedColor)
            {
                vertex.LastColor = vertex.Color;
                vertex.Color = VertexHighlightedColor;
            }

            MainWindow.HighlightedVertices.Add(vertex);
        }

        public static void HighlightEdge(string source, string destination, bool defaultColors, Color customColor)
        {
                ///My.Highlight
                Color VertexHighlightedColor;
                Color EdgeHighlightedColor;
                if (defaultColors)
                {
                    VertexHighlightedColor = (Color)Application.Current.Resources["VertexHighlighted"];
                    EdgeHighlightedColor = (Color)Application.Current.Resources["EdgeHighlighted"];
                }
                else
                {
                    ///VertexHighlightedColor = EdgeHighlightedColor = GetSequentialColor(1, false);
                    VertexHighlightedColor = EdgeHighlightedColor = customColor;
                }

                var vertex1 = GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == source) as TaggedVertex;
                if (vertex1.Color != VertexHighlightedColor)
                {
                    vertex1.LastColor = vertex1.Color;
                    vertex1.Color = VertexHighlightedColor;
                }

                var vertex2 = GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == destination) as TaggedVertex;
                if (vertex2.Color != VertexHighlightedColor)
                {
                    vertex2.LastColor = vertex2.Color;
                    vertex2.Color = EdgeHighlightedColor;
                }

                var edge = GraphToVisualize.Edges.ToList().Find(x => x.Source.ToString() == source && x.Target.ToString() == destination) as TaggedEdge;
                if (edge.Color != EdgeHighlightedColor)
                {
                    edge.LastColor = edge.Color;
                    edge.Color = EdgeHighlightedColor;
                }

                MainWindow.HighlightedEdges.Add(edge);
            

                ///Graph#.Highlight
                /*object vertex1 = GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == source);
                if (!graphLayout.IsHighlightedVertex(vertex1))
                    graphLayout.HighlightVertex(vertex1, null);
                object vertex2 = graphLayout.Graph.Vertices.ToList().Find(x => x.ToString() == destination);
                if (!graphLayout.IsHighlightedVertex(vertex2))
                    graphLayout.HighlightVertex(vertex2, null);
                WeightedEdge<object> edge = graphLayout.Graph.Edges.ToList().Find(x => x.Source.ToString() == source && x.Target.ToString() == destination);
                if (!graphLayout.IsHighlightedEdge(edge))
                    graphLayout.HighlightEdge(edge, null);*/
        }

        public static void HighlightPath(ref int[] path, int pathType, int totalPaths, int timeout, bool cumulative, bool defaultColors, Color customColor)    // ,ref Tuple<int,int>[] ER
            {
                /// Path types: 
                ///     0: the path is [0 or 1, ...] and references to an ER elements (e.g. cycles)
                ///     1: the path is [n € N, ...] and references to ER indexes 
                ///     2: the path is [v1,v2,v3] and references Vertices (e.g. cliques)

                Color VertexHighlightedColor;
                Color EdgeHighlightedColor;
                if (defaultColors)
                {
                    VertexHighlightedColor = (Color)Application.Current.Resources["VertexHighlighted"];
                    EdgeHighlightedColor = (Color)Application.Current.Resources["EdgeHighlighted"];
                }
                else
                {
                    ///VertexHighlightedColor = EdgeHighlightedColor = GetSequentialColor(1, false);
                    VertexHighlightedColor = EdgeHighlightedColor = customColor;
                }

                bool customHighlight = true;
                if (customHighlight)
                {
                   if(!cumulative) MainWindow.RemoveCustomHighlights();
                }
                /*else
                    DataGrid1.IsHitTestVisible = false;*/

                if (timeout == 0)
                {
                    /// Determine the timeout
                    if(totalPaths < 20)  timeout = 6000/totalPaths;
                    else if(totalPaths < 40)  timeout = 8000/totalPaths;
                        else if(totalPaths < 80)  timeout = 12000/totalPaths;
                            else if(totalPaths < 120) timeout = 16000/totalPaths;
                                else if(totalPaths > 120) timeout = 20000/totalPaths;
                }
                else
                {
                    /// Convert from seconds to milliseconds
                    timeout *= 1000; 
                }

                /// For each element in the path
                for(int i = 0; i < path.Length; ++i)
                {
                    if ( ( pathType == 0 && path[i] == 1) || (pathType == 1) )
                    {
                        int idx = (pathType==0)? i: path[i];
                        string v1 = GraphFromLibrary.EdgesReference[idx].Item1.ToString();
                        string v2 = GraphFromLibrary.EdgesReference[idx].Item2.ToString();
                        HighlightEdge(v1, v2, defaultColors, EdgeHighlightedColor);      
                    }

                    if ( pathType == 2 )
                    {
                        if (path.Length > 1)
                        {
                            int idx = path[i];
                            string v1 = idx.ToString();

                            foreach (var adjList in GraphFromLibrary.adjList[idx])
                            {
                                string v2 = adjList.Item1.ToString();

                                /// If v2 is adjacent to v1 AND the v1 is less than v2 (prevents re-processing v2,v1)
                                if ((path.Contains(Convert.ToInt32(v2))) && (Convert.ToInt32(v1) < Convert.ToInt32(v2)))
                                {
                                    HighlightEdge(v1, v2, defaultColors, EdgeHighlightedColor);
                                }
                            }
                        }
                        else
                        {
                            HighlightVertex(path[i].ToString(), defaultColors, VertexHighlightedColor);      
                        }
                    }
                }
                
                ///Grraph#.Highlight
                /*if(!customHighlight)
                {
                    /// Force GUI Refresh
                    MainWindow.ForceRefresh(graphLayout);
                
                    try
                    {
                        /// Display until timeout
                        Thread.Sleep(timeout);

                        /// Or, display while mouse is pressed down
                        //while (MouseDownOnGrid)
                        //{
                            //Console.WriteLine("MouseDownOnGrid is " + MouseDownOnGrid);
                        //}
                    }
                    catch (ThreadInterruptedException ex)
                    {
                        Console.WriteLine("Thread interrupted! (GraphLayout) - " + ex);
                    }

                    /// Restore
                    foreach (var vertex in GraphToVisualize.Vertices)
                        graphLayout.RemoveHighlightFromVertex(vertex);

                    foreach (var edge in GraphToVisualize.Edges)
                        graphLayout.RemoveHighlightFromEdge(edge);
                    DataGrid1.IsHitTestVisible = true;
                }*/
            }

        /// Highlight simplex as colored polygon
        public static void HighlightSimplex(PointCollection myPointCollection, Color color)
        {
            if(myPointCollection.Count == 3)
            {
                Polygon myPolygon = new Polygon();
                myPolygon.Name = "simplex_" + myPointCollection.GetHashCode().ToString();
                myPolygon.Points = myPointCollection;
                myPolygon.Opacity = 0.20;
                myPolygon.Fill = new SolidColorBrush(color);
                myPolygon.Stroke = Brushes.Transparent;
                myPolygon.StrokeThickness = 2;
                //myPolygon1.Stretch = Stretch.Fill;
                var VisualOffset = ((graphLayout.GetType()).GetProperty("VisualOffset", BindingFlags.Instance | BindingFlags.NonPublic)).GetValue(graphLayout, null);
                (myPolygon.GetType()).GetProperty("VisualOffset", BindingFlags.Instance | BindingFlags.NonPublic).SetValue(myPolygon, VisualOffset);
            
                //graphLayout.Children.Add(myPolygon);
                graphLayout.Children.Insert(0, myPolygon);
            }
            else if (myPointCollection.Count > 3)
            {
                ///HighlightSimplex(myPointCollection, color);
            }
        }

        /// Remove highlighted simplices
        public static void RemoveSimplices()
        {
            /// Remove filtrations
            for (int index = graphLayout.Children.Count - 1; index >= 0; index--)
            {
                var children = graphLayout.Children[index];
                if (children is Polygon)
                {
                    Polygon simplex = children as Polygon;
                    if(simplex.Name.StartsWith("simplex_"))
                        graphLayout.Children.RemoveAt(index);
                }
            }
        }

        /// Animated step-by-step filtration
        public static double AnimatedFiltration(int timeSec, int untilFilter)
        {
            /// TODO: When adding/removing, the graphLayout stays unchanged (vertical rectangle), while the graph itself gets smaller
            /// TODO: Canvas size is defined by _topLeft and _bottomRight vertices

            double lastRank = 0;

            Color vertexDisabledColor = (Color) Application.Current.Resources["GraphBackground"];
            Color vertexEnabledColor = HexToColor("#FFFFFF");

            int nfilters = GraphFromLibrary.PH_Filtration.Count;
            if (nfilters > 0)
            {        
                int delayMs = (untilFilter == GraphFromLibrary.PH_Filtration.Count)? (timeSec * 1000 / nfilters) : (0);
                int count = 0;

                /// Draw out all vertices and edges
                foreach (TaggedVertex vertex in GraphToVisualize.Vertices.ToList())
                {
                    ///TODO: Change highlight colors, for vertices and edges
                    
                    vertex.LastColor = vertex.Color;
                    vertex.Color = vertexDisabledColor;

                    /// Hide Vertices' label and circle
                    VertexControl VC = graphLayout.GetVertexControl(vertex);
                    VC.Foreground =  new SolidColorBrush(vertexDisabledColor);
                    VC.BorderBrush =  new SolidColorBrush(vertexDisabledColor);
                    
                    /// Remove highlighted simplices
                    RemoveSimplices();
                }
                foreach (TaggedEdge edge in GraphToVisualize.Edges.ToList())
                {
                    edge.LastColor = edge.Color;
                    edge.Color = vertexDisabledColor;
                }

                ///Refresh & Delay
                if(delayMs>0)
                {
                    ForceRefresh(graphLayout);
                    Delay(1000);
                }

                ///Add filtrations
                if(untilFilter > 0)
                {
                    foreach (var rankedCliques in GraphFromLibrary.PH_Filtration)
                    {
                        if(count < untilFilter)
                        {
                            var filter = rankedCliques.Key;
                            var cliqueS = rankedCliques.Value;
                            lastRank = filter;

                            /// Set color for current filtration (little "darker")
                            Color _color_ = GetSequentialColor(count++, false, 0);
                            var _color_R = (_color_.R * 0.8 > 0)?_color_.R * 0.8 : _color_.R;
                            var _color_G = (_color_.G * 0.8 > 0)?_color_.G * 0.8 : _color_.G;
                            var _color_B = (_color_.B * 0.8 > 0)?_color_.B * 0.8 : _color_.B;
                            Color _color = Color.FromArgb(_color_.A, (byte)_color_R, (byte)_color_G, (byte)_color_B);

                            /// Draw out all previous filtrations
                            foreach (TaggedVertex vertex in GraphToVisualize.Vertices.ToList())
                            {
                                if (vertex.Color != (Color) Application.Current.Resources["GraphBackground"])
                                {
                                    vertex.Color = (Color) Application.Current.Resources["VertexNormal"];

                                    ///Unhide
                                    VertexControl VC = graphLayout.GetVertexControl(vertex);
                                    VC.Foreground =  new SolidColorBrush(vertexEnabledColor);
                                    VC.BorderBrush =  new SolidColorBrush(vertexEnabledColor);
                                }
                            }
                            foreach (TaggedEdge edge in GraphToVisualize.Edges.ToList())
                            {
                                if(edge.Color != (Color)Application.Current.Resources["GraphBackground"])
                                    edge.Color = (Color)Application.Current.Resources["EdgeNormal"];
                            }

                            /// Draw current filtration
                            foreach (var clique in cliqueS)
                            {
                                /// Add simplex
                                var simplex = clique.Item1;
                                PointCollection myPointCollection = new PointCollection();

                                /// Generate "lighter" colors
                                var color_R = (_color.R * 1.25 < 255)?_color.R * 1.25 : _color.R;
                                var color_G = (_color.G * 1.25 < 255)?_color.G * 1.25 : _color.G;
                                var color_B = (_color.B * 1.25 < 255)?_color.B * 1.25 : _color.B;
                                Color color = Color.FromArgb(_color.A, (byte)color_R, (byte)color_G, (byte)color_B);

                                foreach (int v1 in simplex)
                                {
                                    ///var vertex = GraphToVisualize.HiddenVertices.ToList().Find(x => x.ToString() == source.ToString()) as TaggedVertex;
                                    var _source = GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == v1.ToString());
                                    TaggedVertex source = _source as TaggedVertex;
                                    if(source != null)
                                    {
                                        ///GraphToVisualize.UnhideVertex(vertex);
                                        source.Color = color;

                                        /// Add to simplex
                                        myPointCollection.Add(graphLayout.VertexPoint(_source, false));
                                    }

                                    foreach (int v2 in simplex)
                                    {
                                        if(v1 != v2)
                                        {
                                            var _destination = GraphToVisualize.Vertices.ToList().Find(x => x.ToString() == v2.ToString());
                                            TaggedVertex destination = _destination as TaggedVertex;
                                            if(destination != null)
                                            {
                                                destination.Color = color;

                                                ///var edge = GraphToVisualize.HiddenEdges.ToList().Find(x => x.Source.ToString() == source.ToString()) as WeightedEdge<object>;
                                                TaggedEdge edge = GraphToVisualize.Edges.ToList().Find(e => (e.Source.ToString() == source.ToString() && e.Target.ToString() == destination.ToString()) ||
                                                                                                     (e.Source.ToString() == destination.ToString() && e.Target.ToString() == source.ToString())) as TaggedEdge;
                                                if(edge != null)
                                                {
                                                    ///GraphToVisualize.UnhideEdge(edge);
                                                    edge.Color = color;
                                                }
                                            }
                                        }
                                    }
                                }
                                if (clique.Item1.Length == 3)
                                    HighlightSimplex(myPointCollection, color);
                            }

                            ///Delay
                            if(delayMs > 0)
                            {
                                ForceRefresh(graphLayout);
                                Delay(delayMs);
                            }
                        }
                    }
                }
            }

            /// Finally, set standard colors for all
            if(untilFilter == GraphFromLibrary.PH_Filtration.Count)
            {
                foreach (TaggedVertex vertex in GraphToVisualize.Vertices.ToList())
                {
                    vertex.Color = (Color)Application.Current.Resources["VertexNormal"];

                    /// Unhide Vertices' label and circle
                    VertexControl VC = graphLayout.GetVertexControl(vertex);
                    VC.Foreground =  new SolidColorBrush(vertexEnabledColor);
                    VC.BorderBrush =  new SolidColorBrush(vertexEnabledColor);
                }
                foreach (TaggedEdge edge in GraphToVisualize.Edges.ToList())
                {
                    edge.Color = (Color)Application.Current.Resources["EdgeNormal"];
                }
            }

            return lastRank;
        }

        #endregion

        #region EVENT HANDLERS

        private void RadioBSource1_Checked(object sender, RoutedEventArgs e)
        {
            /// File import
            ShowHideExtraFields(false);

            ComboGraph.Items.Clear();
            //string inputDirectory = Application.Current.MainWindow.Resources["InputDirectory"].ToString();    //not yet available during startup
            string inputDirectory = System.IO.Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location) + @"\input\";
            var here = Directory.EnumerateFiles(inputDirectory, "input_*.txt");
            int itemSelected = 0;
            foreach (string file in here)
            {
                string fileName = System.IO.Path.GetFileName(file);
                ComboGraph.Items.Add(fileName);
                if (fileName == "input_sample_0008_stella.txt")
                    itemSelected = ComboGraph.Items.Count-1;
            }

            ComboGraph.SelectedIndex = itemSelected;
        }

        private void RadioBSource2_Checked(object sender, RoutedEventArgs e)
        {
            /// Hypercube
            ShowHideExtraFields(false);

            ComboGraph.Items.Clear();
            ComboGraph.Items.Add("Hypercube n:4 m: 4");
            ComboGraph.Items.Add("Hypercube n:8 m: 12");
            ComboGraph.Items.Add("Hypercube n:16 m: 32");
            ComboGraph.Items.Add("Hypercube n:32 m: 80");
            ComboGraph.Items.Add("Hypercube n:64 m: 192");
            ComboGraph.Items.Add("Hypercube n:128 m: 448");
            ComboGraph.Items.Add("Hypercube n:256 m: 1024");
            ComboGraph.Items.Add("Hypercube n:512 m: 2304");
            ComboGraph.SelectedIndex = 0;
        }

        private void RadioBSource4_Checked(object sender, RoutedEventArgs e)
        {
            /// Density graph
            ShowHideExtraFields(true);
        }

        private void RadioBSource3_Checked(object sender, RoutedEventArgs e)
        {
            /// Euclidean graph
            ShowHideExtraFields(true);
        }

        private void RadioBAlgo1_Checked(object sender, RoutedEventArgs e)
        {
            ComboAlgorithm.Items.Clear();
            ComboAlgorithm.Items.Add("Horton");
            ComboAlgorithm.Items.Add("De Pina");
            ComboAlgorithm.Items.Add("Hybrid");
            ComboAlgorithm.Items.Add("Amaldi");
            ComboAlgorithm.SelectedIndex = 3;

            /// Load available parameters
            SetParamsOnGUI();
        }

        private void RadioBAlgo2_Checked(object sender, RoutedEventArgs e)
        {
            ComboAlgorithm.Items.Clear();
            ComboAlgorithm.Items.Add("BronKerbosch (naive)");
            ComboAlgorithm.Items.Add("BronKerbosch (tomita)");
            ComboAlgorithm.Items.Add("BronKerbosch (eppstein)");
            ComboAlgorithm.Items.Add("Criccaldi");
            ComboAlgorithm.SelectedIndex = 2;

            /// Load available parameters
            SetParamsOnGUI();
        }

        private void RadioBAlgo3_Checked(object sender, RoutedEventArgs e)
        {
            ComboAlgorithm.Items.Clear();

            /*for (int i = 8; i <= 10; ++i)
            {
                ComboAlgorithm.Items.Add(GraphFromLibrary.algorithmS[i]);
            }*/
            
            ComboAlgorithm.Items.Add("Javaplex");
            ComboAlgorithm.Items.Add("Gudhi");
            ComboAlgorithm.Items.Add("PHAT");
            ///ComboAlgorithm.Items.Add("DiPHA");
            ///ComboAlgorithm.Items.Add("Perseus");
            /// 
            ComboAlgorithm.SelectedIndex = 0;

            /// Load available parameters
            SetParamsOnGUI();
        }

        private void ComboGraph_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            /// Hypercube graph
            ShowHideExtraFields(false);

            /// Prevent building graph during "Items.clear" event
            if(ComboGraph.SelectedItem != null)
                /// Set graph creation as "pending"
                pendingChanges = true;
        }

        private void SliderA_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            if(SliderA != null && SliderB != null )
            {
                int V = Convert.ToInt32(SliderA.Value);
                double D = SliderB.Value;
                int maxE = V * (V - 1) / 2;
                int E = (int)(maxE * D);

                /// euclidean graph
                if (RadioBSource3.IsChecked == true)    
                {  
                    SliderB.Maximum = maxE;
                }

                /// density graph
                if (RadioBSource4.IsChecked == true)    
                {  
                    double temp = SliderB.Value;
                    SliderB.Value = 0;
                    SliderB.Value = temp;
                    Label_Bval.Content = "(" + (E).ToString() + " edges)";
                }
            }

            /// Set graph creation as "pending"
            pendingChanges = true;
        }

        private void SliderB_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
                int V = Convert.ToInt32(SliderA.Value);
                double D = e.NewValue;
                int maxE = V * (V - 1) / 2;
                int E = (int)(maxE * D);
            
                ///density graph
                if (RadioBSource4.IsChecked == true)
                {
                    /// update #edges
                    Label_Bval.Content = "(" + (E).ToString() + " edges)";
                }

                /// Set graph creation as "pending"
                pendingChanges = true;
        }

        private void author_click(object sender, MouseEventArgs e)
        {
            System.Diagnostics.Process.Start("https://github.com/mitxael");
        }
        
        private void ComboAlgorithm_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            /// Load available parameters
            SetParamsOnGUI();
        }

        #endregion

        #region BGWORKER       
        /// CALLER:
        //backgroundWorker1.RunWorkerAsync();

        private void InitializeBackgroundWorker()
        {
            backgroundWorker1.DoWork += new DoWorkEventHandler(backgroundWorker1_DoWork);
            backgroundWorker1.RunWorkerCompleted += new RunWorkerCompletedEventHandler(backgroundWorker1_RunWorkerCompleted);
            backgroundWorker1.ProgressChanged += new ProgressChangedEventHandler(backgroundWorker1_ProgressChanged);
            backgroundWorker1.WorkerReportsProgress = true;
        }

        private void backgroundWorker1_DoWork(object sender, DoWorkEventArgs e)
        {
            /// RUN ALL BACKGROUND TASKS HERE
            
            /// Get the BackgroundWorker that raised this event.
            BackgroundWorker worker = sender as BackgroundWorker;

            /// Assign the result of the computation to the Result property of the DoWorkEventArgs object
            //e.Result = ComputeFibonacci((int)e.Argument, worker, e);
        }

        private void backgroundWorker1_RunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
        {
            ///UPDATE UI WHEN THE WORK IS COMPLETE
            
            if (e.Error != null)
            {
                /// First, handle the case where an exception was thrown.
                MessageBox.Show(e.Error.Message);
            }
            else if (e.Cancelled)
            {
                /// Next, handle the case where the user canceled the operation.
                //resultLabel.Text = "Canceled";
            }
            else
            {
                /// Finally, handle the case where the operation succeeded.
                //resultLabel.Text = e.Result.ToString();
            }
        }
        
        private void backgroundWorker1_ProgressChanged(object sender, ProgressChangedEventArgs e)
        {
            /// IF NEEDED, UPDATE AFTER CHANGES
            //this.progressBar1.Value = e.ProgressPercentage;
        }

        #endregion

   
    }

    #region OTHER-CLASSES

    [ValueConversion(typeof(object), typeof(string))]
    public class StringConverter : IValueConverter
    {
        public object Convert(object value, Type targetType, object parameter, CultureInfo culture)
        {
            return value == null ? null : value.ToString();
        }

        public object ConvertBack(object value, Type targetType, object parameter, CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }

    public class EdgeColorConverter : IValueConverter
    {

        public object Convert(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
        {
            return new SolidColorBrush((Color) value);
        }

        public object ConvertBack(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
        {
            throw new NotImplementedException();
        }
    }

    public class VertexColorConverter : IValueConverter
    {

        public object Convert(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
        {
            return new SolidColorBrush((Color) value);
        }

        public object ConvertBack(object value, Type targetType, object parameter, System.Globalization.CultureInfo culture)
        {
            throw new NotImplementedException();
        }

    }

    #endregion
}