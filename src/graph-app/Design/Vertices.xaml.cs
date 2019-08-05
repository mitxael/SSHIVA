using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using CHoleR.Helper;
using GraphSharp;
using GraphSharp.Controls;
using GraphSharp.Controls.Zoom;


namespace CHoleR.Design
{
    partial class Vertices : ResourceDictionary
    { 
        public Vertices()
        {
            InitializeComponent();
        }    
        
        public void ChangeVertexColor_OnClick(object sender, RoutedEventArgs e)
        {
            Color highlightColor = (Color) Application.Current.Resources["VertexHighlighted"];

            var vertexcontrol = sender as VertexControl;
            var vertex = vertexcontrol.DataContext as TaggedVertex;
            var graphlayout = vertexcontrol.Parent as GraphLayout;
            var zoomcontrol = graphlayout.Parent as ZoomControl;
            
            if (vertex is TaggedVertex)
            {
               if(vertex.Color != highlightColor)
               {
                   vertex.LastColor = vertex.Color;
                   vertex.Color = highlightColor;

                   ///Change also its edges
                   /*var edges = graphlayout.Graph.Edges.Where( x =>
                       ((TaggedVertex)x.Source == vertex) || ((TaggedVertex) x.Target == vertex));
                   foreach (TaggedEdge edge in edges)
                   {
                       edge.Color = Colors.Red;
                   }*/
               }
               else
               {
                   vertex.Color = vertex.LastColor;
                   vertex.LastColor = highlightColor;

                   ///Change also its edges
                   /*var edges = graphlayout.Graph.Edges.Where( x =>
                       ((TaggedVertex)x.Source == vertex) || ((TaggedVertex) x.Target == vertex));
                   foreach (TaggedEdge edge in edges)
                   {
                       edge.Color = Colors.Red;
                   }*/
               }
            }
        }

        public void ContextualMenu(object sender, RoutedEventArgs e)
        {
            MenuItem menuitem = sender as MenuItem;
            Vertex vertex = menuitem.DataContext as Vertex;

            switch (menuitem.Tag.ToString())
            {
                case "1":
                    HighlightCycles_OnClick(vertex);
                    break;
                case "2":
                    HighlightCliques_OnClick(vertex);
                    break;
                case "3":
                    HideVertex_OnClick(vertex);                
                    break;
            }
        }

        public void HideVertex_OnClick(Vertex vertex)
        {           
            if (vertex is Vertex)
            {
                ///Remove vertex
                //MainWindow.GraphToVisualize.RemoveVertex(vertex);
                MainWindow.GraphToVisualize.HideVertex(vertex);
            }
        }
        
        public void UnhideVertex_OnClick(Vertex vertex)
        {           
            if (vertex is Vertex)
            {               
                ///Unhide ALL (add vertex and saved edges)
                MainWindow.GraphToVisualize.UnhideVertexAndEdges(vertex);
            }
        }

        public void HighlightCycles_OnClick(Vertex vertex)
        {            
            if (vertex is Vertex)
            {
                ///Highlight cycles
                MainWindow.RemoveCustomHighlights();

                int v = Convert.ToInt32(vertex.ToString());
                int count = 0;

                foreach (var cycle in MainWindow.GraphFromLibrary.minimumCycleBasis)
                {
                    foreach (var edgeRef in cycle.Select((value, idx) => new { value, idx }))
                    {
                        var edge = MainWindow.GraphFromLibrary.EdgesReference[edgeRef.idx];
                        if ( (edgeRef.value == 1) && (edge.Item1 == v || edge.Item2 == v) )
                        {
                            int[] data = cycle;
                            Color color = MainWindow.GetSequentialColor(count++, false, 0);
                            MainWindow.HighlightPath(ref data, 0, 1, 1, true, false, color);
                            break;
                        }
                    }
                }

                /// Highlight selected vertex
                //Color vcolor = (Color)Application.Current.Resources["VertexHighlighted"];
                Color vcolor = Colors.Black;
                MainWindow.HighlightVertex(vertex.ToString(), false, vcolor);
            }
        }

        public void HighlightCliques_OnClick(Vertex vertex)
        {
            if (vertex is Vertex)
            {
                ///Highlight ciques
                MainWindow.RemoveCustomHighlights();
                 
                int count = 0;
                foreach (var clique in MainWindow.GraphFromLibrary.maximalCliqueSet)
                {
                    if(clique.ToList().Contains(Convert.ToInt32(vertex.ToString())))
                    {
                        int[] data = clique;
                        Color color = MainWindow.GetSequentialColor(count++, false, 0);
                        MainWindow.HighlightPath(ref data, 2, 1, 1, true, false, color);
                    }
                }

                /// Highlight selected vertex
                //Color vcolor = (Color)Application.Current.Resources["VertexHighlighted"];
                Color vcolor = Colors.Black;
                MainWindow.HighlightVertex(vertex.ToString(), false, vcolor);
            }
        }        
    }
}