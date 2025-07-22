import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import 'package:image_picker/image_picker.dart';
import 'dart:io';
import '../services/payment_service.dart';
import '../services/analysis_service.dart';

class PhotoPicker extends StatefulWidget {
  @override
  _PhotoPickerState createState() => _PhotoPickerState();
}

class _PhotoPickerState extends State<PhotoPicker> {
  final _urlController = TextEditingController();
  final _imagePicker = ImagePicker();
  
  @override
  Widget build(BuildContext context) {
    return Consumer2<PaymentService, AnalysisService>(
      builder: (context, paymentService, analysisService, child) {
        return Container(
          width: double.infinity,
          padding: EdgeInsets.all(20),
          decoration: BoxDecoration(
            color: Colors.white.withOpacity(0.08),
            borderRadius: BorderRadius.circular(8),
          ),
          // FLATTENED STRUCTURE - Row instead of Column to prevent overflow
          child: Row(
            children: [
              // Upload button - compact
              Expanded(
                flex: 1,
                child: _buildUploadButton(paymentService, analysisService),
              ),
              
              SizedBox(width: 20),
              
              // URL input section - compact
              Expanded(
                flex: 2,
                child: _buildUrlSection(paymentService, analysisService),
              ),
            ],
          ),
        );
      },
    );
  }
  
  Widget _buildUploadButton(PaymentService paymentService, AnalysisService analysisService) {
    return GestureDetector(
      onTap: () => _pickImage(paymentService, analysisService),
      child: Container(
        height: 100, // Fixed height prevents overflow
        decoration: BoxDecoration(
          border: Border.all(
            color: Colors.white.withOpacity(0.7),
            style: BorderStyle.solid,
            width: 1,
          ),
          borderRadius: BorderRadius.circular(8),
        ),
        child: Column(
          mainAxisAlignment: MainAxisAlignment.center,
          mainAxisSize: MainAxisSize.min, // KEY: Prevents overflow
          children: [
            CustomPaint(
              size: Size(24, 24),
              painter: UploadIconPainter(),
            ),
            SizedBox(height: 8),
            Text(
              'Upload',
              style: TextStyle(
                color: Colors.white.withOpacity(0.7),
                fontSize: 14,
                fontWeight: FontWeight.w500,
              ),
            ),
          ],
        ),
      ),
    );
  }
  
  Widget _buildUrlSection(PaymentService paymentService, AnalysisService analysisService) {
    return Container(
      height: 100, // Fixed height prevents overflow
      child: Column(
        mainAxisAlignment: MainAxisAlignment.center,
        children: [
          Text(
            'Or paste image URL:',
            style: TextStyle(
              color: Colors.white.withOpacity(0.7),
              fontSize: 12,
            ),
          ),
          SizedBox(height: 8),
          Row(
            children: [
              Expanded(
                child: TextField(
                  controller: _urlController,
                  decoration: InputDecoration(
                    hintText: 'Image URL...',
                    isDense: true,
                    contentPadding: EdgeInsets.symmetric(horizontal: 12, vertical: 8),
                  ),
                  style: TextStyle(color: Colors.white, fontSize: 14),
                  onSubmitted: (_) => _analyzeUrl(paymentService, analysisService),
                ),
              ),
              SizedBox(width: 8),
              GestureDetector(
                onTap: () => _analyzeUrl(paymentService, analysisService),
                child: Container(
                  padding: EdgeInsets.all(8),
                  decoration: BoxDecoration(
                    color: Colors.white.withOpacity(0.08),
                    borderRadius: BorderRadius.circular(8),
                  ),
                  child: CustomPaint(
                    size: Size(16, 16),
                    painter: SearchIconPainter(),
                  ),
                ),
              ),
            ],
          ),
        ],
      ),
    );
  }
  
  Future<void> _pickImage(PaymentService paymentService, AnalysisService analysisService) async {
    if (!paymentService.canAnalyze()) {
      paymentService.showPayment();
      return;
    }
    
    try {
      final XFile? image = await _imagePicker.pickImage(
        source: ImageSource.gallery,
        maxWidth: 1920,
        maxHeight: 1080,
        imageQuality: 85,
      );
      
      if (image != null) {
        final imageFile = File(image.path);
        await analysisService.analyzeImage(imageFile);
        await paymentService.incrementUsage();
        
        ScaffoldMessenger.of(context).showSnackBar(
          SnackBar(content: Text('Image analysis complete')),
        );
      }
    } catch (error) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('Image pickup failed: $error')),
      );
    }
  }
  
  Future<void> _analyzeUrl(PaymentService paymentService, AnalysisService analysisService) async {
    final url = _urlController.text.trim();
    if (url.isEmpty) return;
    
    if (!paymentService.canAnalyze()) {
      paymentService.showPayment();
      return;
    }
    
    try {
      _urlController.clear();
      await analysisService.analyzeImageUrl(url);
      await paymentService.incrementUsage();
      
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('URL analysis complete')),
      );
    } catch (error) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(content: Text('URL analysis failed: $error')),
      );
    }
  }
  
  @override
  void dispose() {
    _urlController.dispose();
    super.dispose();
  }
}

// Custom painter for upload icon matching original SVG
class UploadIconPainter extends CustomPainter {
  @override
  void paint(Canvas canvas, Size size) {
    final paint = Paint()
      ..color = Colors.white.withOpacity(0.7)
      ..style = PaintingStyle.stroke
      ..strokeWidth = 2.0;

    final scale = size.width / 24.0;
    
    // Cloud base: path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"
    final cloudPath = Path();
    cloudPath.moveTo(21 * scale, 15 * scale);
    cloudPath.lineTo(21 * scale, 19 * scale);
    cloudPath.quadraticBezierTo(21 * scale, 21 * scale, 19 * scale, 21 * scale);
    cloudPath.lineTo(5 * scale, 21 * scale);
    cloudPath.quadraticBezierTo(3 * scale, 21 * scale, 3 * scale, 19 * scale);
    cloudPath.lineTo(3 * scale, 15 * scale);
    
    canvas.drawPath(cloudPath, paint);
    
    // Upload arrow: polyline points="7,10 12,15 17,10"
    final arrowPath = Path();
    arrowPath.moveTo(7 * scale, 10 * scale);
    arrowPath.lineTo(12 * scale, 15 * scale);
    arrowPath.lineTo(17 * scale, 10 * scale);
    
    canvas.drawPath(arrowPath, paint);
    
    // Upload line: line x1="12" y1="15" x2="12" y2="3"
    canvas.drawLine(
      Offset(12 * scale, 15 * scale),
      Offset(12 * scale, 3 * scale),
      paint,
    );
  }

  @override
  bool shouldRepaint(CustomPainter oldDelegate) => false;
}

// Custom painter for search icon matching original SVG
class SearchIconPainter extends CustomPainter {
  @override
  void paint(Canvas canvas, Size size) {
    final paint = Paint()
      ..color = Colors.white
      ..style = PaintingStyle.stroke
      ..strokeWidth = 2.0;

    final scale = size.width / 24.0;
    
    // Magnifying glass: circle cx="11" cy="11" r="8"
    canvas.drawCircle(
      Offset(11 * scale, 11 * scale),
      8 * scale,
      paint,
    );
    
    // Handle: path d="M21 21l-4.35-4.35"
    canvas.drawLine(
      Offset(21 * scale, 21 * scale),
      Offset(16.65 * scale, 16.65 * scale),
      paint,
    );
  }

  @override
  bool shouldRepaint(CustomPainter oldDelegate) => false;
} 