import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import 'package:flutter_svg/flutter_svg.dart';
import 'payment_service.dart';
import 'analysis_service.dart';

class InputBar extends StatefulWidget {
  @override
  _InputBarState createState() => _InputBarState();
}

class _InputBarState extends State<InputBar> {
  final _controller = TextEditingController();
  
  @override
  Widget build(BuildContext context) {
    return Consumer2<PaymentService, AnalysisService>(
      builder: (context, paymentService, analysisService, child) {
        return Row(
          children: [
            Expanded(
              child: TextField(
                controller: _controller,
                decoration: InputDecoration(
                  hintText: 'Molecule name or formula...',
                  border: OutlineInputBorder(
                    borderRadius: BorderRadius.circular(8),
                    borderSide: BorderSide(color: Colors.white.withOpacity(0.3)),
                  ),
                  enabledBorder: OutlineInputBorder(
                    borderRadius: BorderRadius.circular(8),
                    borderSide: BorderSide(color: Colors.white.withOpacity(0.3)),
                  ),
                  focusedBorder: OutlineInputBorder(
                    borderRadius: BorderRadius.circular(8),
                    borderSide: BorderSide(color: Colors.white),
                  ),
                  filled: true,
                  fillColor: Colors.white.withOpacity(0.08),
                  hintStyle: TextStyle(color: Colors.white.withOpacity(0.7)),
                  contentPadding: EdgeInsets.symmetric(horizontal: 16, vertical: 12),
                ),
                style: TextStyle(color: Colors.white, fontSize: 16),
                onSubmitted: (value) => _analyzeText(paymentService, analysisService, value),
              ),
            ),
            SizedBox(width: 12),
            GestureDetector(
              onTap: () {
                if (!paymentService.canAnalyze()) {
                  paymentService.showPayment();
                }
              },
              child: Container(
                padding: EdgeInsets.all(12),
                decoration: BoxDecoration(
                  color: Colors.white.withOpacity(0.08),
                  borderRadius: BorderRadius.circular(8),
                ),
                child: SvgPicture.asset(
                  'assets/account.svg',
                  width: 20,
                  height: 20,
                  colorFilter: ColorFilter.mode(
                    Colors.white.withOpacity(0.7),
                    BlendMode.srcIn,
                  ),
                ),
              ),
            ),
          ],
        );
      },
    );
  }
  
  Future<void> _analyzeText(PaymentService paymentService, AnalysisService analysisService, String text) async {
    if (text.trim().isEmpty) return;
    
    if (!paymentService.canAnalyze()) {
      paymentService.showPayment();
      return;
    }
    
    try {
      _controller.clear();
      await analysisService.analyzeText(text.trim());
      await paymentService.incrementUsage();
      
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Analysis complete')),
      );
    } catch (error) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Analysis failed: $error')),
      );
    }
  }
  
  @override
  void dispose() {
    _controller.dispose();
    super.dispose();
  }
} 