using System;
using System.Windows.Data;
using System.Diagnostics;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Windows;
using System.Windows.Media;
using Brush = System.Windows.Media.Brush;
using Point = System.Windows.Point;
using Size = System.Windows.Size;

namespace GraphSharp.Converters
{
	/// <summary>
	/// Converts the position and sizes of the source and target points, and the route informations
	/// of an edge to a path.
	/// The edge can bend, or it can be straight line.
	/// </summary>
	public class EdgeRouteToPathConverter : IMultiValueConverter
	{
		public object Convert( object[] values, Type targetType, object parameter, System.Globalization.CultureInfo culture )
		{
			Debug.Assert( values != null && values.Length == 10, "EdgeRouteToPathConverter should have 10 parameters: pos (1,2), size (3,4) of source; pos (5,6), size (7,8) of target; routeInformation (9); showArrows (10)." );

			#region Get the inputs
			//get the position of the source
			var sourcePos = new Point
								{
									X = ( values[0] != DependencyProperty.UnsetValue ? (double)values[0] : 0.0 ),
									Y = ( values[1] != DependencyProperty.UnsetValue ? (double)values[1] : 0.0 )
								};
			//get the size of the source
			var sourceSize = new Size
								{
									Width = ( values[2] != DependencyProperty.UnsetValue ? (double)values[2] : 0.0 ),
									Height = ( values[3] != DependencyProperty.UnsetValue ? (double)values[3] : 0.0 )
								};
			//get the position of the target
			var targetPos = new Point
								{
									X = ( values[4] != DependencyProperty.UnsetValue ? (double)values[4] : 0.0 ),
									Y = ( values[5] != DependencyProperty.UnsetValue ? (double)values[5] : 0.0 )
								};
			//get the size of the target
			var targetSize = new Size
								{
									Width = ( values[6] != DependencyProperty.UnsetValue ? (double)values[6] : 0.0 ),
									Height = ( values[7] != DependencyProperty.UnsetValue ? (double)values[7] : 0.0 )
								};

			//get the route informations
			Point[] routeInformation = ( values[8] != DependencyProperty.UnsetValue ? (Point[])values[8] : null );
			//get showArrows
			Boolean showArrows = (values[9] != DependencyProperty.UnsetValue ? (Boolean)values[9] : true);
			#endregion

			bool hasRouteInfo = routeInformation != null && routeInformation.Length > 0;


			/// Create the path
			Point p1 = GraphConverterHelper.CalculateAttachPoint( sourcePos, sourceSize, ( hasRouteInfo ? routeInformation[0] : targetPos ) );
			/// Create the arrows
			Point p2 = GraphConverterHelper.CalculateAttachPoint( targetPos, targetSize, ( hasRouteInfo ? routeInformation[routeInformation.Length - 1] : sourcePos ) );

			var segments = new PathSegment[1 + ( hasRouteInfo ? routeInformation.Length : 0 )];
			if ( hasRouteInfo )
			{
				///append route points
				for ( int i = 0; i < routeInformation.Length; i++ )
					segments[i] = new LineSegment( routeInformation[i], true );
			}
			Point pLast = ( hasRouteInfo ? routeInformation[routeInformation.Length - 1] : p1 );
			Vector v = pLast - p2;
			v = v / v.Length * 5;
			Vector n = new Vector( -v.Y, v.X ) * 0.3;

			if (showArrows == true)
			{ 
				segments[segments.Length - 1] = new LineSegment( p2 + v, true );

				return new PathFigureCollection(2)
				{
				new PathFigure(p1, segments, false),
				
				new PathFigure(p2,
							   new PathSegment[]
							   {
								   new LineSegment(p2 + v - n, true),
								   new LineSegment(p2 + v + n, true)
							   }, true)
				};
			}
			else
			{
				segments[segments.Length - 1] = new LineSegment( p2 , true );

				return new PathFigureCollection(2)
				{
					new PathFigure(p1, segments, false),

					/// Add weight. Pass the "middle point" of the edge
					//EdgeWeight(new Point( p1.X + 0.5*(p2.X-p1.X), p1.Y + 0.5*(p2.Y-p1.Y)), n)
				};
			}
		}

		public object[] ConvertBack( object value, Type[] targetTypes, object parameter, System.Globalization.CultureInfo culture )
		{
			throw new NotSupportedException();
		}


		/// *************** UTILITIES ***********************************************************

		public PathFigure EdgeWeight(Point startP, Vector dev)
		{
		    return new PathFigure(startP, new PathSegment[]
			{
				new LineSegment(startP - dev, true)
			}, true);
		}
	}
}