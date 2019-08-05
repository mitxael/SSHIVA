using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Controls;
using System.Windows.Documents;
using System.Windows.Media;

namespace CHoleR.Helper
{
    public static class Logger
    {
        public static RichTextBox LogBox { get; set; }

        public static void WriteLogInfo(string info)
        {
            TextRange tr = new TextRange(LogBox.Document.ContentEnd, LogBox.Document.ContentEnd);
            string TimePattern = DateTime.Now.ToString("dd/MM/yy HH:mm:ss.fff") + " | >>> ";

            tr.Text = TimePattern + info + "\n";
            tr.ApplyPropertyValue(TextElement.ForegroundProperty, Brushes.Black);
        }

        public static void WriteLogException(string message)
        {
            TextRange tr = new TextRange(LogBox.Document.ContentStart, LogBox.Document.ContentStart);
            string TimePattern = DateTime.Now.ToString("dd/MM/yy HH:mm:ss.fff") + " | >>> ";

            tr.Text = TimePattern + message + "\n";
            tr.ApplyPropertyValue(TextElement.ForegroundProperty, Brushes.Red);
        }

        public static void ClearLog()
        {
            LogBox.Document.Blocks.Clear();
        }
    }
}
