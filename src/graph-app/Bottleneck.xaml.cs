using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using GRAPH;
using Microsoft.Win32;

namespace CHoleR
{
    /// <summary>
    /// Interaction logic for Bottleneck.xaml
    /// </summary>
    public partial class Bottleneck : Window
    {
        private string file1 = "";
        private string file2 = "";
        /// LastPath's can be a List<string> based on selectedindex
        private string lastPath1 = Application.Current.MainWindow.Resources[$"InputDirectory"].ToString();
        private string lastPath2 = Application.Current.MainWindow.Resources[$"InputDirectory"].ToString();
        private bool InitializationCompleted = false;

        public Bottleneck()
        {
            InitializeComponent();
        }

        #region BUTTON ACTIONS

        /// Compute Bottleneck Distance
        private async void BtnComputeBottlebeck_Click(object sender, RoutedEventArgs e)
        { 
            double b = 0;
            /**string inputDirectory = Application.Current.MainWindow.Resources[$"InputDirectory"].ToString();

            string file1 = inputDirectory + "persistance_diag1.txt";
            string file2 = inputDirectory + "persistance_diag2.txt";**/

            /// EXECUTE FUNCTION
            HomologyManager ph = new HomologyManager();
            try
            {
                ///Run task    
                Task t = Task.Run(() => b = ph.bottleneckDistance(file1, file2, 0));
                await Task.WhenAll(t);
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.Message);
                throw;
            }
            finally
            {
                if (b == -1)
                {
                    MessageBox.Show("ERROR. File1 is invalid!");
                }
                else if (b == -2)
                {
                    MessageBox.Show("ERROR. File2 is invalid!");
                }
                else
                {
                    Label_Result.Content = b;
                    Label_Result.Visibility = Visibility.Visible;    
                }
            }
        }

        private void BtnHelp_Click(object sender, RoutedEventArgs e)
        {
            MessageBox.Show("Input files must contain in each line the \n" +
                            "birth and death of classes. For instance:\n" +
                            "\n" +
                            "\t2 7\n" +
                            "\t5 11\n" +
                            "\t8 inf\n" +
                            "\n" +
                            "Optionally, each line can also specify the value\n" +
                            "of the error bound.\n");
        }

        #endregion

        #region EVENTS

        /// Startup actions
        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            /// Set Window Position
            var desktopWorkingArea = System.Windows.SystemParameters.WorkArea;
            this.Left = desktopWorkingArea.Left;
            this.Top = desktopWorkingArea.Bottom - this.Height;
            
            /// Set Window Background
            this.Background = new SolidColorBrush( (Color)Application.Current.Resources["GraphBackground"] );

            /// Set combo #1
            var dimensions = MainWindow.GraphFromLibrary.PH_barcodes.Keys.ToArray();
            foreach (int dim in dimensions)
            {
                ComboPD1.Items.Add("Graph PH from dimension-" + dim);
            }
            ComboPD1.Items.Add("<browse for a file...>");
            ComboPD1.SelectedIndex = 0;

            /// Set combo #2
            ComboPD2.Items.Add("...");
            ComboPD2.Items.Add("Browse for a file...");
            ComboPD2.SelectedIndex = 0;

            /// Everything was successfully initialized
            InitializationCompleted = true;
        }

        /// Combo1 selection changed
        private void ComboPD1_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            int lastItem = ComboPD1.Items.Count - 1;
            if (ComboPD1.SelectedIndex == lastItem)
            {
                OpenFileDialog openFileDialog = new OpenFileDialog();
                openFileDialog.InitialDirectory = Application.Current.MainWindow.Resources[$"InputDirectory"].ToString();
                if(openFileDialog.ShowDialog() == true)
                {
                    ComboPD1.Items.Insert(lastItem, openFileDialog.SafeFileName); ///File.ReadAllText(openFileDialog.FileName));
                    file1 = openFileDialog.FileName;
                    lastPath1 = openFileDialog.InitialDirectory;
                    ComboPD1.SelectedIndex = ComboPD1.Items.Count - 2;
                }
            }
            else
            {
                var dimensions = MainWindow.GraphFromLibrary.PH_barcodes.Keys.ToArray();
                if (ComboPD1.SelectedIndex < dimensions.Length)
                {
                    file1 = System.IO.Path.GetTempPath() + "persistence.txt";   ///C:\Users\UserName\AppData\Local\Temp\persistence.txt
                    string content1 = "";
                    var intervals = MainWindow.GraphFromLibrary.PH_barcodes[ComboPD1.SelectedIndex];
                    foreach (var pair in intervals)
                    {
                        string birth = pair.Key.Item1.ToString();
                        string death = (pair.Key.Item2 < Double.PositiveInfinity)? pair.Key.Item2.ToString() : "inf";
                        content1 += birth + " " + death + Environment.NewLine;
                    }
                    File.WriteAllText(file1, content1);
                }
                else
                {
                    file1 = lastPath1 + ComboPD1.Items[ComboPD1.SelectedIndex];
                }
            }
        }
        
        /// Combo2 selection changed
        private void ComboPD2_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            
            if(InitializationCompleted)
            {
                int lastItem = ComboPD2.Items.Count - 1;
                if(ComboPD2.SelectedIndex == lastItem)
                {
                    OpenFileDialog openFileDialog = new OpenFileDialog();
                    openFileDialog.InitialDirectory = Application.Current.MainWindow.Resources[$"InputDirectory"].ToString();
                    if (openFileDialog.ShowDialog() == true)
                    {
                        ComboPD2.Items.Insert(lastItem, openFileDialog.SafeFileName); ///File.ReadAllText(openFileDialog.FileName));
                        file2 = openFileDialog.FileName;
                        lastPath2 = openFileDialog.InitialDirectory;
                        ComboPD2.SelectedIndex =lastItem;
                    }
                }
                else
                {
                    file2 = lastPath2 + ComboPD2.Items[ComboPD2.SelectedIndex];
                }
            }
        }

        #endregion
    }
}
