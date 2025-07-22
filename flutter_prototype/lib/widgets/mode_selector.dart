import 'package:flutter/material.dart';
import 'package:provider/provider.dart';
import '../services/payment_service.dart';
import '../screens/home_screen.dart'; // For AnalysisMode enum

class ModeSelector extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Consumer<PaymentService>(
      builder: (context, paymentService, child) {
        return Container(
          padding: EdgeInsets.all(8),
          decoration: BoxDecoration(
            color: Colors.white.withOpacity(0.08), // Interactive element background
            borderRadius: BorderRadius.circular(8),
          ),
          child: Row(
            children: [
              Expanded(
                child: _buildModeTab(
                  context,
                  paymentService,
                  AnalysisMode.camera,
                  _buildCameraIcon(),
                  'Camera',
                ),
              ),
              SizedBox(width: 8),
              Expanded(
                child: _buildModeTab(
                  context,
                  paymentService,
                  AnalysisMode.photo,
                  _buildPhotoIcon(),
                  'Photo',
                ),
              ),
            ],
          ),
        );
      },
    );
  }

  // Original camera icon: simple lens (two concentric circles)
  Widget _buildCameraIcon() {
    return CustomPaint(
      size: Size(16, 16),
      painter: CameraLensIconPainter(),
    );
  }

  // Original photo SVG from HTML: landscape/mountain scene
  Widget _buildPhotoIcon() {
    return CustomPaint(
      size: Size(16, 16),
      painter: PhotoIconPainter(),
    );
  }
  
  Widget _buildModeTab(
    BuildContext context,
    PaymentService paymentService,
    AnalysisMode mode,
    Widget icon,
    String label,
  ) {
    final isSelected = paymentService.selectedMode == mode;
    
    return GestureDetector(
      onTap: () => paymentService.setMode(mode),
      child: Container(
        padding: EdgeInsets.symmetric(vertical: 12, horizontal: 16),
        decoration: BoxDecoration(
          color: isSelected 
            ? Colors.white.withOpacity(0.08) 
            : Colors.transparent,
          borderRadius: BorderRadius.circular(8),
        ),
        child: Row(
          mainAxisAlignment: MainAxisAlignment.center,
          children: [
            Theme(
              data: Theme.of(context).copyWith(
                iconTheme: IconThemeData(
                  color: isSelected 
                    ? Colors.white 
                    : Colors.white.withOpacity(0.7),
                ),
              ),
              child: icon,
            ),
            SizedBox(width: 8),
            Text(
              label,
              style: TextStyle(
                color: isSelected 
                  ? Colors.white 
                  : Colors.white.withOpacity(0.7),
                fontSize: 14,
                fontWeight: FontWeight.w500,
              ),
            ),
          ],
        ),
      ),
    );
  }
}

// Simple camera lens icon - matches original lens.svg
class CameraLensIconPainter extends CustomPainter {
  @override
  void paint(Canvas canvas, Size size) {
    final paint = Paint()
      ..color = Colors.white
      ..style = PaintingStyle.stroke
      ..strokeWidth = 2.0;

    final fillPaint = Paint()
      ..color = Colors.white.withOpacity(0.8)
      ..style = PaintingStyle.fill;

    final scale = size.width / 24.0;
    final center = Offset(12 * scale, 12 * scale);
    
    // Inner lens (filled circle): cx="12" cy="12" r="5" fill="#ffffff" opacity="0.8"
    canvas.drawCircle(
      center,
      5 * scale,
      fillPaint,
    );
    
    // Outer lens ring (stroke circle): cx="12" cy="12" r="9" stroke="#ffffff" fill="none"
    canvas.drawCircle(
      center,
      9 * scale,
      paint,
    );
  }

  @override
  bool shouldRepaint(CustomPainter oldDelegate) => false;
}

// Custom painter for photo icon matching original SVG
class PhotoIconPainter extends CustomPainter {
  @override
  void paint(Canvas canvas, Size size) {
    final paint = Paint()
      ..color = Colors.white
      ..style = PaintingStyle.stroke
      ..strokeWidth = 2.0;

    final scale = size.width / 24.0;
    
    // Photo frame: rect x="3" y="3" width="18" height="18" rx="2" ry="2"
    canvas.drawRRect(
      RRect.fromLTRBR(
        3 * scale, 3 * scale, 21 * scale, 21 * scale,
        Radius.circular(2 * scale)
      ),
      paint,
    );
    
    // Sun/moon: circle cx="9" cy="9" r="2"
    canvas.drawCircle(
      Offset(9 * scale, 9 * scale),
      2 * scale,
      paint,
    );
    
    // Mountain/landscape path: M21 15l-3.086-3.086a2 2 0 0 0-2.828 0L6 21
    final mountainPath = Path();
    mountainPath.moveTo(21 * scale, 15 * scale);
    mountainPath.lineTo(17.914 * scale, 11.914 * scale);
    mountainPath.quadraticBezierTo(
      16.5 * scale, 10.5 * scale,
      15.086 * scale, 11.914 * scale
    );
    mountainPath.lineTo(6 * scale, 21 * scale);
    
    canvas.drawPath(mountainPath, paint);
  }

  @override
  bool shouldRepaint(CustomPainter oldDelegate) => false;
} 