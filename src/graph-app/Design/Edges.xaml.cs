using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Media;
using CHoleR.Helper;
using GraphSharp;
using GraphSharp.Controls;
using GraphSharp.Controls.Zoom;


namespace CHoleR.Design
{
    partial class Edges : ResourceDictionary
    { 
        public Edges()
        {
            InitializeComponent();
        }

        public void ChangeEdgeColor_OnClick(object sender, RoutedEventArgs e)
        {
            Color highlightColor = (Color) Application.Current.Resources["EdgeHighlighted"];

            var edgecontrol = sender as EdgeControl;
            var edge = edgecontrol.DataContext as TaggedEdge;
            var graphlayout = edgecontrol.Parent as GraphLayout;
            var zoomcontrol = graphlayout.Parent as ZoomControl;

            if (edge is TaggedEdge)
            {
                if (edge.Color != highlightColor)
                {
                    edge.LastColor = edge.Color;
                    edge.Color = highlightColor;
                }
                else
                {
                    edge.Color = edge.LastColor;
                    edge.LastColor = highlightColor;
                }
                
            }
        }
    }
}