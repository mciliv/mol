import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import '../services/payment_service.dart';
import '../services/analysis_service.dart';

class InputBar extends StatefulWidget {
  @override
  _InputBarState createState() => _InputBarState();
}

class _InputBarState extends State<InputBar> {
  final _textController = TextEditingController();
  
  @override
  Widget build(BuildContext context) {
    return Consumer2<PaymentService, AnalysisService>(
      builder: (context, paymentService, analysisService, child) {
        return Container(
          padding: EdgeInsets.all(16),
          decoration: BoxDecoration(
            color: Colors.white.withOpacity(0.08), // Interactive element background
            borderRadius: BorderRadius.circular(8),
          ),
          child: Row(
            children: [
              // Text input
              Expanded(
                child: TextField(
                  controller: _textController,
                  decoration: InputDecoration(
                    hintText: 'Type any molecule name...',
                    border: InputBorder.none,
                    isDense: true,
                    contentPadding: EdgeInsets.symmetric(horizontal: 16, vertical: 12),
                  ),
                  style: TextStyle(color: Colors.white),
                  onSubmitted: (text) => _handleSubmit(context, paymentService, analysisService),
                ),
              ),
              
              SizedBox(width: 16),
              
              // Account button with original credit card SVG
              _buildAccountButton(paymentService),
            ],
          ),
        );
      },
    );
  }
  
  Widget _buildAccountButton(PaymentService paymentService) {
    return GestureDetector(
      onTap: () => _handleAccountTap(paymentService),
      child: Container(
        padding: EdgeInsets.symmetric(horizontal: 12, vertical: 8),
        decoration: BoxDecoration(
          color: Colors.white.withOpacity(0.08),
          borderRadius: BorderRadius.circular(8),
        ),
        child: Row(
          mainAxisSize: MainAxisSize.min,
          children: [
            // Original credit card SVG icon
            CustomPaint(
              size: Size(16, 16),
              painter: CreditCardIconPainter(
                color: _getAccountButtonColor(paymentService),
              ),
            ),
            SizedBox(width: 8),
            Text(
              _getAccountButtonText(paymentService),
              style: TextStyle(
                color: _getAccountButtonColor(paymentService),
                fontWeight: FontWeight.w500,
                fontSize: 14,
              ),
            ),
          ],
        ),
      ),
    );
  }
  
  Color _getAccountButtonColor(PaymentService paymentService) {
    if (paymentService.hasPaymentMethod) {
      return Colors.green; // Payment setup
    } else {
      return Colors.orange; // Needs setup
    }
  }
  
  String _getAccountButtonText(PaymentService paymentService) {
    if (paymentService.hasPaymentMethod) {
      return paymentService.userName ?? 'Account';
    } else {
      return 'Add Card';
    }
  }
  
  void _handleAccountTap(PaymentService paymentService) {
    if (!paymentService.hasPaymentMethod) {
      paymentService.showPayment();
    } else {
      // Show account management (placeholder)
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(
          content: Text('Usage: ${paymentService.usageCount} analyses'),
          duration: Duration(seconds: 2),
        ),
      );
    }
  }
  
  Future<void> _handleSubmit(
    BuildContext context,
    PaymentService paymentService,
    AnalysisService analysisService,
  ) async {
    final text = _textController.text.trim();
    if (text.isEmpty) return;
    
    // Check payment method
    if (!paymentService.canAnalyze()) {
      paymentService.showPayment();
      return;
    }
    
    // Clear input
    _textController.clear();
    
    // Analyze
    try {
      await analysisService.analyzeText(text);
      await paymentService.incrementUsage();
      
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Analysis complete for: $text')),
      );
    } catch (error) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Analysis failed: $error')),
      );
    }
  }
  
  @override
  void dispose() {
    _textController.dispose();
    super.dispose();
  }
}

// Custom painter for credit card icon matching original SVG
class CreditCardIconPainter extends CustomPainter {
  final Color color;
  
  CreditCardIconPainter({required this.color});
  
  @override
  void paint(Canvas canvas, Size size) {
    final paint = Paint()
      ..color = color
      ..style = PaintingStyle.stroke
      ..strokeWidth = 2.0;

    final scale = size.width / 24.0;
    
    // Credit card: rect x="2" y="4" width="20" height="16" rx="2"
    canvas.drawRRect(
      RRect.fromLTRBR(
        2 * scale, 4 * scale, 22 * scale, 20 * scale,
        Radius.circular(2 * scale)
      ),
      paint,
    );
    
    // Card stripe: line x1="2" y1="10" x2="22" y2="10"
    canvas.drawLine(
      Offset(2 * scale, 10 * scale),
      Offset(22 * scale, 10 * scale),
      paint,
    );
  }

  @override
  bool shouldRepaint(CreditCardIconPainter oldDelegate) => color != oldDelegate.color;
} 